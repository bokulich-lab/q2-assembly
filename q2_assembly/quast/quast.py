# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import json
import os
import subprocess
import tempfile
from distutils.dir_util import copy_tree
from typing import List, Union

import pandas as pd
import pkg_resources
import q2templates
from q2_types.per_sample_sequences import \
    (SingleLanePerSamplePairedEndFastqDirFmt,
     SingleLanePerSampleSingleEndFastqDirFmt)
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt

from .._utils import (run_command, _remove_html_element, _modify_links,
                      _construct_param, _process_common_input_params)

TEMPLATES = pkg_resources.resource_filename('q2_assembly', 'assets')


def _get_sample_from_path(fp):
    """Extracts sample name from a contig's file path."""
    return os.path.basename(fp).split('_', maxsplit=1)[0]


def _process_quast_arg(arg_key, arg_val):
    """Creates a list with argument and its value.

    Argument values represented by a list will be converted to a single
    string joined by commas, e.g.: [1, 2, 3] -> '1,2,3'.
    Argument names will be converted to command line parameters by
    appending a '--' prefix and replacing all '_' with '-',
    e.g.: 'some_parameter' -> '--some-parameter'.

    Args:
        arg_key (str): Argument name.
        arg_val: Argument value.

    Returns:
        (converted_arg, arg_value): Tuple containing a prepared command line
            parameter and its value.
    """
    if arg_key == 'threads' and (not arg_val or arg_val > 1):
        # TODO: this needs to be fixed (to allow multiprocessing)
        print('Multiprocessing is currently not supported. Resetting '
              'number of threads to 1.')
        arg_value = '1'
    elif not isinstance(arg_val, list):
        arg_value = str(arg_val)
    else:
        arg_value = ','.join(str(x) for x in arg_val)
    return _construct_param(arg_key), arg_value


def _evaluate_contigs(
        results_dir: str, contigs: ContigSequencesDirFmt, reads: dict,
        paired: bool, common_args: list
) -> List[str]:
    # TODO: this will probably get replaced by "quast.py" once we
    #  want to extend the support beyond metagenomes
    cmd = ['metaquast.py', '-o', results_dir]
    cmd.extend(common_args)
    samples = []

    for fp in sorted(glob.glob(os.path.join(str(contigs), '*_contigs.fa'))):
        cmd.append(fp)
        samples.append(_get_sample_from_path(fp))

    if reads:
        rev_count = sum(
            [True if x['rev'] else False for _, x in reads.items()]
        )
        if (rev_count < len(samples) > rev_count) and paired:
            raise Exception(
                f'Number of reverse reads ({rev_count}) does not match '
                f'the number of provided contig files ({len(samples)}). '
                f'Please check your input files.')

        for s in samples:
            try:
                if paired:
                    cmd.extend(
                        ['--pe1', reads[s].get('fwd'),
                         '--pe2', reads[s].get('rev')]
                    )
                else:
                    cmd.extend(
                        ['--single', reads[s].get('fwd')]
                    )
            except KeyError:
                # TODO: improve this msg by adding name of the missing sample
                raise Exception(
                    'Some samples are missing from the reads file. '
                    'Please check your input files.'
                )

    try:
        run_command(cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception('An error was encountered while running '
                        f'QUAST, (return code {e.returncode}), '
                        'please inspect stdout and stderr to learn more.')

    return samples


def _fix_html_reports(results_dp: str):
    # remove link to icarus browser from report.html
    # (otherwise, the visualisation can get messed up once clicked)
    report_fp = os.path.join(results_dp, 'report.html')
    _remove_html_element(report_fp, 'p', elem_id='icarus')

    # remove "Main menu" button from contig browser
    # (otherwise, the visualisation can get messed up once clicked)
    contig_browser_fp = os.path.join(
        results_dp, 'icarus_viewers', 'contig_size_viewer.html')
    _remove_html_element(
        contig_browser_fp, 'div', elem_id='to_main_menu_button')

    # make all the external links open in a new tab
    _modify_links(report_fp)


def evaluate_contigs(
        # TODO: expose more parameters
        output_dir: str,
        contigs: ContigSequencesDirFmt,
        reads: Union[SingleLanePerSamplePairedEndFastqDirFmt,
                     SingleLanePerSampleSingleEndFastqDirFmt] = None,
        min_contig: int = None,
        threads: int = None,
        k_mer_stats: bool = False,
        k_mer_size: int = None,
        contig_thresholds: List[int] = None,
        x_for_Nx: int = None
):

    common_args = _process_common_input_params(
        processing_func=_process_quast_arg,
        min_contig=min_contig, threads=threads, k_mer_stats=k_mer_stats,
        k_mer_size=k_mer_size, contig_thresholds=contig_thresholds,
        x_for_Nx=x_for_Nx
    )

    reads_fps = {}
    paired = False
    if reads:
        paired = isinstance(reads, SingleLanePerSamplePairedEndFastqDirFmt)
        manifest = reads.manifest.view(pd.DataFrame)
        for samp in list(manifest.index):
            reads_fps[samp] = {
                'fwd': manifest.loc[samp, 'forward'],
                'rev': manifest.loc[samp, 'reverse'] if paired else None
            }

    with tempfile.TemporaryDirectory() as tmp:
        results_dir = os.path.join(tmp, 'results')

        # run quast
        samples = _evaluate_contigs(
            results_dir, contigs, reads_fps, paired, common_args)

        # fix/remove some URLs
        _fix_html_reports(results_dir)

        copy_tree(os.path.join(TEMPLATES, 'quast'), output_dir)
        copy_tree(results_dir, os.path.join(output_dir, 'quast_data'))

        context = {
            'tabs': [
                {
                    'title': 'QC report',
                    'url': 'index.html'
                },
                {
                    'title': 'Contig browser',
                    'url': 'q2_icarus.html'

                },
            ],
            'samples': json.dumps(samples)
        }

        index = os.path.join(TEMPLATES, 'quast', 'index.html')
        icarus = os.path.join(TEMPLATES, 'quast', 'q2_icarus.html')

        templates = [index, icarus]
        q2templates.render(templates, output_dir, context=context)
