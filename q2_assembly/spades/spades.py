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

import pandas as pd
from q2_types.per_sample_sequences import (
    ContigSequencesDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)

from .._utils import (
    _construct_param,
    _process_common_input_params,
    concatenate_files,
    get_file_extension,
    run_command,
)


def _process_spades_arg(arg_key, arg_val):
    """Creates a list with argument and its value.

    Argument values represented by a list will be converted to a single
    string joined by spaces, e.g.: [1, 2, 3] -> '1 2 3'.
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
    elif arg_key == "phred_offset" and arg_val == "auto-detect":
        return []
    elif not isinstance(arg_val, list):
        return [_construct_param(arg_key), str(arg_val)]
    elif arg_key == "k":
        return ["-k", ",".join(str(x) for x in arg_val)]
    else:
        arg_value = ",".join(str(x) for x in arg_val)
        return [_construct_param(arg_key), arg_value]


def _process_sample(sample, fwd, rev, common_args, out):
    """Constructs assembly command for SPAdes and runs the assembly.

    Args:
        sample: Name of the sample to be processed.
        fwd: Location of the forward reads.
        rev: Location of the reverse reads. If set to None, single-end reads
            are assumed.
        common_args: list of arguments that should be passed to the spades
            assembly command
        out: output format
    """
    with tempfile.TemporaryDirectory() as tmp:
        results_dir = os.path.join(tmp, "results")
        cmd = ["spades.py"]
        if rev:
            cmd.extend(["-1", fwd, "-2", rev])
        else:
            cmd.extend(["-s", fwd])
        cmd.extend(["-o", results_dir])
        cmd.extend(common_args)

        try:
            run_command(cmd, verbose=True)
        except subprocess.CalledProcessError as e:
            raise Exception(
                "An error was encountered while running "
                f"SPAdes, (return code {e.returncode}), "
                "please inspect stdout and stderr to learn more."
            )

        shutil.move(
            os.path.join(results_dir, "contigs.fasta"),
            os.path.join(str(out), f"{sample}_contigs.fa"),
        )


def _assemble_spades(
    seqs, meta, common_args, coassemble=False
) -> ContigSequencesDirFmt:
    """Runs the assembly for all available samples.

    Both, paired- and single-end reads can be processed - the output will
    be adjusted accordingly.

    Args:
        seqs: Sequences to be processed.
        meta: True if it is a metagenomic assembly.
        common_args: List of common flags and their values for
            the metaSPAdes command.
        coassemble: True if user wants to coassemble reads
            from all samples.

    Returns:
        result (ContigSequencesDirFmt): Assembled contigs.

    """

    paired = isinstance(seqs, SingleLanePerSamplePairedEndFastqDirFmt)
    if not paired and meta:
        raise NotImplementedError(
            'SPAdes v3.15.2 in "meta" mode supports only ' "paired-end reads."
        )
    manifest = seqs.manifest.view(pd.DataFrame)
    result = ContigSequencesDirFmt()

    if coassemble:
        print(
            "WARNING: The co-assembly process includes concatenating contents of all "
            "read files. This process may take a while to finish!"
        )

        fwds = manifest["forward"]  # this will be  a list
        revs = (
            manifest["reverse"] if paired else None
        )  # this will be a list if it exists

        with tempfile.TemporaryDirectory() as tmpdir:
            extension = get_file_extension(fwds[0])
            fwd = os.path.join(tmpdir, f"all_contigs_fwd{extension}")
            rev = (
                os.path.join(tmpdir, f"all_contigs_rev{extension}") if paired else None
            )

            concatenate_files(fwds, fwd)
            if paired:
                concatenate_files(revs, rev)

            _process_sample("all_contigs", fwd, rev, common_args, result)

    else:
        for samp in list(manifest.index):
            fwd = manifest.loc[samp, "forward"]
            rev = manifest.loc[samp, "reverse"] if paired else None

            _process_sample(samp, fwd, rev, common_args, result)
    return result


def assemble_spades(
    seqs: Union[
        SingleLanePerSamplePairedEndFastqDirFmt, SingleLanePerSampleSingleEndFastqDirFmt
    ],
    isolate: bool = False,
    sc: bool = False,
    meta: bool = False,
    bio: bool = False,
    corona: bool = False,
    plasmid: bool = False,
    metaviral: bool = False,
    metaplasmid: bool = False,
    only_assembler: bool = False,
    careful: bool = False,
    disable_rr: bool = False,
    threads: int = 1,
    memory: int = 250,
    k: List[int] = ["auto"],
    cov_cutoff: Union[float, str] = "off",
    phred_offset: str = "auto-detect",
    debug: bool = False,
    coassemble: bool = False,
) -> ContigSequencesDirFmt:
    kwargs = {k: v for k, v in locals().items() if k not in ["seqs", "coassemble"]}
    common_args = _process_common_input_params(
        processing_func=_process_spades_arg, params=kwargs
    )

    return _assemble_spades(
        seqs=seqs, meta=meta, coassemble=coassemble, common_args=common_args
    )
