# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import shutil
import subprocess
import tempfile
from typing import List, Union

import pandas as pd
from q2_types.per_sample_sequences import \
    (SingleLanePerSamplePairedEndFastqDirFmt,
     SingleLanePerSampleSingleEndFastqDirFmt)
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt

from .._utils import (run_command, _construct_param,
                      _process_common_input_params)


def _process_megahit_arg(arg_key, arg_val):
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
    if not isinstance(arg_val, list):
        arg_value = str(arg_val)
    else:
        arg_value = ','.join(str(x) for x in arg_val)
    return _construct_param(arg_key), arg_value


def _process_sample(sample, fwd, rev, common_args, out):
    """Constructs assembly command for MEGAHIT and runs the assembly.

    Args:
        sample: Name of the sample to be processed.
        fwd: Location of the forward reads.
        rev: Location of the reverse reads. If set to None, single-end reads
            are assumed.
        common_args: list of arguments that should be passed to the megahit
            assembly command
        out: output format
    """
    with tempfile.TemporaryDirectory() as tmp:
        results_dir = os.path.join(tmp, 'results')
        cmd = ['megahit']
        if rev:
            cmd.extend(['-1', fwd, '-2', rev])
        else:
            cmd.extend(['-r', fwd])
        cmd.extend(['-o', results_dir])
        cmd.extend(common_args)

        try:
            run_command(cmd, verbose=True)
        except subprocess.CalledProcessError as e:
            raise Exception('An error was encountered while running MEGAHIT, '
                            f'(return code {e.returncode}), please inspect '
                            'stdout and stderr to learn more.')

        shutil.move(os.path.join(results_dir, 'final.contigs.fa'),
                    os.path.join(str(out), f'{sample}_contigs.fa'))


def _assemble_megahit(
        seqs, presets, min_count, k_list, k_min, k_max, k_step, no_mercy,
        bubble_level, prune_level, prune_depth, disconnect_ratio,
        low_local_ratio, max_tip_len, cleaning_rounds, no_local, kmin_1pass,
        memory, mem_flag, num_cpu_threads, no_hw_accel, min_contig_len
) -> ContigSequencesDirFmt:
    """Runs the assembly for all available samples.

    Both, paired- and single-end reads can be processed - the output will
    be adjusted accordingly.

    Returns:
        result (ContigSequencesDirFmt): Assembled contigs.

    """
    common_args = _process_common_input_params(
        processing_func=_process_megahit_arg,
        presets=presets, min_count=min_count, k_list=k_list,
        k_min=k_min, k_max=k_max, k_step=k_step, no_mercy=no_mercy,
        bubble_level=bubble_level, prune_level=prune_level,
        prune_depth=prune_depth, disconnect_ratio=disconnect_ratio,
        low_local_ratio=low_local_ratio, max_tip_len=max_tip_len,
        cleaning_rounds=cleaning_rounds, no_local=no_local,
        kmin_1pass=kmin_1pass, memory=memory, mem_flag=mem_flag,
        num_cpu_threads=num_cpu_threads, no_hw_accel=no_hw_accel,
        min_contig_len=min_contig_len
    )

    paired = isinstance(seqs, SingleLanePerSamplePairedEndFastqDirFmt)
    manifest = seqs.manifest.view(pd.DataFrame)
    result = ContigSequencesDirFmt()

    for samp in list(manifest.index):
        fwd = manifest.loc[samp, 'forward']
        rev = manifest.loc[samp, 'reverse'] if paired else None

        _process_sample(samp, fwd, rev, common_args, result)
    return result


def assemble_megahit(
        seqs: Union[SingleLanePerSamplePairedEndFastqDirFmt,
                    SingleLanePerSampleSingleEndFastqDirFmt],
        presets: str = None,
        min_count: int = None,
        k_list: List[int] = None,
        k_min: int = None,
        k_max: int = None,
        k_step: int = None,
        no_mercy: bool = False,
        bubble_level: int = None,
        prune_level: int = None,
        prune_depth: int = None,
        disconnect_ratio: float = None,
        low_local_ratio: float = None,
        max_tip_len: int = None,
        cleaning_rounds: int = None,
        no_local: bool = False,
        kmin_1pass: bool = False,
        memory: float = None,
        mem_flag: int = None,
        num_cpu_threads: int = None,
        no_hw_accel: bool = False,
        min_contig_len: int = None
) -> ContigSequencesDirFmt:

    return _assemble_megahit(
        seqs=seqs, presets=presets, min_count=min_count, k_list=k_list,
        k_min=k_min, k_max=k_max, k_step=k_step, no_mercy=no_mercy,
        bubble_level=bubble_level, prune_level=prune_level,
        prune_depth=prune_depth, disconnect_ratio=disconnect_ratio,
        low_local_ratio=low_local_ratio, max_tip_len=max_tip_len,
        cleaning_rounds=cleaning_rounds, no_local=no_local,
        kmin_1pass=kmin_1pass, memory=memory, mem_flag=mem_flag,
        num_cpu_threads=num_cpu_threads, no_hw_accel=no_hw_accel,
        min_contig_len=min_contig_len)
