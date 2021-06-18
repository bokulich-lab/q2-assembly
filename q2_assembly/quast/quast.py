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

from .._utils import run_command, _remove_html_element, _modify_links

TEMPLATES = pkg_resources.resource_filename('q2_assembly', 'assets')


def _get_sample_from_path(fp):
    """Extracts sample name from a contig's file path."""
    return os.path.basename(fp).split('_', maxsplit=1)[0]


def _evaluate_contigs(
        results_dir: str, contigs: ContigSequencesDirFmt, reads: dict,
        paired: bool, min_contig: int, threads: int
) -> List[str]:
    # TODO: this will probably get replaced by "quast.py" once we
    #  want to extend the support beyond metagenomes
    cmd = ['metaquast.py', '-o', results_dir]
    samples = []

    if min_contig:
        cmd.extend(['-m', str(min_contig)])

    if not threads or threads > 1:
        # TODO: this needs to be fixed (to allow multiprocessing)
        print('Multiprocessing is currently not supported. Resetting '
              'number of threads to 1.')
    cmd.extend(['-t', str(1)])

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
        threads: int = None
):
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
            results_dir, contigs, reads_fps, paired, min_contig, threads)

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
