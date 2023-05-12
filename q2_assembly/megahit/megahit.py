# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
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
import warnings

import pandas as pd
from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt

from .._utils import _construct_param, _process_common_input_params, run_command


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
        [converted_arg, arg_value]: List containing a prepared command line
            parameter and its value.
    """
    if isinstance(arg_val, bool) and arg_val:
        return [_construct_param(arg_key)]
    elif not isinstance(arg_val, list):
        return [_construct_param(arg_key), str(arg_val)]
    else:
        arg_value = ",".join(str(x) for x in arg_val)
        return [_construct_param(arg_key), arg_value]


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
        results_dir = os.path.join(tmp, "results")
        cmd = ["megahit"]
        if rev:
            cmd.extend(["-1", fwd, "-2", rev])
        else:
            cmd.extend(["-r", fwd])
        cmd.extend(["-o", results_dir])
        cmd.extend(common_args)

        try:
            run_command(cmd, verbose=True)
        except subprocess.CalledProcessError as e:
            raise Exception(
                "An error was encountered while running MEGAHIT, "
                f"(return code {e.returncode}), please inspect "
                "stdout and stderr to learn more."
            )

        shutil.move(
            os.path.join(results_dir, "final.contigs.fa"),
            os.path.join(str(out), f"{sample}_contigs.fa"),
        )


def _assemble_megahit(seqs, common_args) -> ContigSequencesDirFmt:
    """Runs the assembly for all available samples.

    Both, paired- and single-end reads can be processed - the output will
    be adjusted accordingly.

    Args:
        seqs: Sequences to be processed.
        common_args: List of common flags and their values for
            the MEGAHIT command.

    Returns:
        result (ContigSequencesDirFmt): Assembled contigs.

    """

    paired = isinstance(seqs, SingleLanePerSamplePairedEndFastqDirFmt)
    manifest = seqs.manifest.view(pd.DataFrame)
    result = ContigSequencesDirFmt()

    for samp in list(manifest.index):
        fwd = manifest.loc[samp, "forward"]
        rev = manifest.loc[samp, "reverse"] if paired else None

        _process_sample(samp, fwd, rev, common_args, result)
    return result


def warn_about_presets():
    warning = (
        "The presets parameter overrides settings for the min_count and k_list "
        "parameters. The settings of min_count and k_list registered in provenance "
        "may not reflect the actual settings used by the presets parameter. Refer "
        "to the megahit documentation for more details, and refer to the presets "
        "values of these parameters when interpreting or reporting your results."
    )
    warnings.warn(warning, UserWarning)


def assemble_megahit(
    seqs: Union[
        SingleLanePerSamplePairedEndFastqDirFmt, SingleLanePerSampleSingleEndFastqDirFmt
    ],
    presets: str = None,
    min_count: int = 2,
    k_list: List[int] = [21, 29, 39, 59, 79, 99, 119, 141],
    k_min: int = None,
    k_max: int = None,
    k_step: int = None,
    no_mercy: bool = False,
    bubble_level: int = 2,
    prune_level: int = 2,
    prune_depth: int = 2,
    disconnect_ratio: float = 0.1,
    low_local_ratio: float = 0.2,
    max_tip_len: int = "auto",
    cleaning_rounds: int = 5,
    no_local: bool = False,
    kmin_1pass: bool = False,
    memory: float = 0.9,
    mem_flag: int = 1,
    num_cpu_threads: int = 1,
    no_hw_accel: bool = False,
    min_contig_len: int = 200,
) -> ContigSequencesDirFmt:
    if max_tip_len == "auto":
        max_tip_len = None
    if presets == "disabled":
        presets = None
    else:
        warn_about_presets()
    if any([k_min, k_max, k_step]) and not all([k_min, k_max, k_step]):
        raise ValueError(
            "If any of the parameters k_min, k_max, or k_step are used "
            "then all must be explicitly set."
        )

    kwargs = {k: v for k, v in locals().items() if k not in ["seqs"]}
    common_args = _process_common_input_params(
        processing_func=_process_megahit_arg, params=kwargs
    )

    return _assemble_megahit(seqs=seqs, common_args=common_args)
