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
from bs4 import BeautifulSoup as BS
from q2_types.per_sample_sequences import \
    (SingleLanePerSamplePairedEndFastqDirFmt,
     SingleLanePerSampleSingleEndFastqDirFmt)
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt

from q2_assembly._utils import run_command

TEMPLATES = pkg_resources.resource_filename('q2_assembly', 'assets')


def _remove_html_element(fp: str, tag: str, elem_id: str):
    """Removes an HTML tag and its contents.

    Uses BeautifulSoup to open an HTML file, find an element by tag and
    its id and remove it, if found. The original file will be overwritten
    by its modified version.

    Args:
         fp (str): Path to the original HTML file.
         tag (str): Type of the tag to be removed.
         elem_id (str): ID of the element (tag) to be removed.
    """
    with open(fp, 'r') as r:
        soup = BS(r.read(), 'html.parser')
        element = soup.find(tag, id=elem_id)
        if element:
            element.decompose()
    os.remove(fp)

    with open(fp, 'w') as r:
        r.write(str(soup))


def _modify_links(fp: str):
    """Modifies all "a" tags to automatically open in a new tab rather
        than the original iFrame.

    Uses BeautifulSoup to open an HTML file, find all "a" tags and
    add a "target" property. The original file will be overwritten
    by its modified version.

    Args:
         fp (str): Path to the original HTML file.
    """
    with open(fp, 'r') as r:
        soup = BS(r.read(), 'html.parser')
        links = soup.find_all('a')
        for line in links:
            line['target'] = '_blank'
    os.remove(fp)

    with open(fp, 'w') as r:
        r.write(str(soup))


def _get_sample_from_path(fp):
    """Extracts sample name from a contig's file path."""
    return os.path.basename(fp).split('_', maxsplit=1)[0]


def _evaluate_contigs(
        results_dir, contigs, reads, paired, min_contig, threads
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

        # remove link to icarus browser from report.html
        # (otherwise, the visualisation can get messed up once clicked)
        report_fp = os.path.join(results_dir, 'report.html')
        _remove_html_element(report_fp, 'p', elem_id='icarus')

        # remove "Main menu" button from contig browser
        # (otherwise, the visualisation can get messed up once clicked)
        contig_browser_fp = os.path.join(
            results_dir, 'icarus_viewers', 'contig_size_viewer.html')
        _remove_html_element(
            contig_browser_fp, 'div', elem_id='to_main_menu_button')

        # make all the external links open in a new tab
        _modify_links(report_fp)

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
