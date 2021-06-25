# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import os
import shutil
import subprocess
import tempfile
from copy import deepcopy
from typing import List

import biom
import pandas as pd
import skbio
from q2_types.feature_data import DNAFASTAFormat
from q2_types.per_sample_sequences import \
    (CasavaOneEightSingleLanePerSampleDirFmt)

from .._utils import (run_command, _process_common_input_params)


def _process_iss_arg(arg_key, arg_val):
    """Creates a list with argument and its value.

    Argument values represented by a list will be converted to a single
    string joined by spaces, e.g.: [1, 2, 3] -> '1 2 3'.
    Argument names will be converted to command line parameters by
    appending a '--' prefix, e.g.: 'some_parameter' -> '--some_parameter'.

    Args:
        arg_key (str): Argument name.
        arg_val: Argument value.

    Returns:
        [converted_arg, arg_value, ...]: List containing a prepared command
            line parameter and its value(s).
    """
    if not isinstance(arg_val, list):
        return [f'--{arg_key}', str(arg_val)]
    else:
        flags = [f'--{arg_key}']
        flags.extend([str(x) for x in arg_val])
        return flags


def _rename_reads_files(dp):
    reads = sorted(glob.glob(os.path.join(dp, '*.fastq.gz')))
    reads_new = [r.replace('.fastq.gz', '_001.fastq.gz') for r in reads]
    for r, rn in zip(reads, reads_new):
        os.rename(r, rn)
    return reads_new


def _generate_reads(samples, args, result_fp):
    base_cmd = ['iss', 'generate', '--compress']
    base_cmd.extend(args)

    for s in samples:
        cmd = deepcopy(base_cmd)
        sample_prefix = os.path.join(result_fp, f'{s}_00_L001')

        cmd.extend(['--output', sample_prefix])

        try:
            run_command(cmd, verbose=True)
        except subprocess.CalledProcessError as e:
            raise Exception(
                'An error was encountered while running InSilicoSeq, '
                f'(return code {e.returncode}), please inspect '
                'stdout and stderr to learn more.')

    _rename_reads_files(result_fp)


def _abundances_to_biom(abundance_fps):
    abundances = []
    for f in abundance_fps:
        sample = os.path.splitext(os.path.basename(f))[0].split('_')[0]
        abundances.append(
            pd.read_csv(f, sep='\t', index_col=0, header=None, names=[sample]))
    biom_df = pd.concat(abundances, axis=1)
    return biom.Table(data=biom_df.values,
                      observation_ids=biom_df.index.tolist(),
                      sample_ids=biom_df.columns.tolist())


# TODO: allow to input custom genomes and sample from those
def generate_reads(
        sample_names: List[str] = None,
        ncbi: List[str] = None,
        n_genomes_ncbi: List[int] = None,
        abundance: str = None,
        coverage: str = None,
        n_reads: int = None,
        mode: str = None,
        model: str = 'HiSeq',
        gc_bias: bool = False,
        cpus: int = None,
        debug: bool = False,
        seed: int = None
) -> (CasavaOneEightSingleLanePerSampleDirFmt, DNAFASTAFormat, biom.Table):

    if n_genomes_ncbi and ncbi and (len(n_genomes_ncbi) != len(ncbi)):
        raise Exception('Genome counts (--n_genomes_ncbi) need to correspond'
                        'to the kingdoms names (--ncbi). You provided '
                        f'{len(ncbi)} kingdom(s) but {len(n_genomes_ncbi)} '
                        f'corresponding genome counts were found. Please '
                        f'correct your input.')

    if len(set(sample_names)) < len(sample_names):
        dupl = set([str(x) for x in sample_names if sample_names.count(x) > 1])
        raise Exception('Sample names need to be unique. Found duplicated '
                        f'names: {", ".join(sorted(dupl))}')

    kwargs = {k: v for k, v in locals().items() if k not in ['sample_names']}
    args = _process_common_input_params(
        processing_func=_process_iss_arg, params=kwargs
    )

    with tempfile.TemporaryDirectory() as tmp:
        result_reads = CasavaOneEightSingleLanePerSampleDirFmt()
        result_genomes = DNAFASTAFormat()

        # simulate reads
        _generate_reads(sample_names, args, tmp)

        # move reads into CasavaFmt
        for f in glob.glob(os.path.join(tmp, '*.fastq.gz')):
            shutil.move(f, str(result_reads))

        # move original genomes into DNAFASTAFmt
        with result_genomes.open() as fout:
            for f in sorted(glob.glob(os.path.join(tmp, '*.fasta'))):
                for seq in skbio.read(f, format='fasta'):
                    seq.write(fout)

        # convert abundances to a biom table
        abund_suffix = 'coverage' if coverage else 'abundance'
        abundance_fps = sorted(
            glob.glob(os.path.join(tmp, f'*_{abund_suffix}.txt')))
        result_biom = _abundances_to_biom(abundance_fps)

        return result_reads, result_genomes, result_biom
