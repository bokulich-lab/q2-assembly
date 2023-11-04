# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
import subprocess
from copy import deepcopy

from q2_types.bowtie2 import Bowtie2IndexDirFmt
from q2_types_genomics.per_sample_data import (
    ContigSequencesDirFmt,
    MultiBowtie2IndexDirFmt,
    MultiMAGSequencesDirFmt,
)

from q2_assembly._utils import _process_common_input_params, run_command
from q2_assembly.bowtie2.utils import (
    _get_subdir_from_path, _process_bowtie2build_arg, _assert_inputs_not_empty
)


def _index_seqs(
    fasta_fps: list, result_fp: str, common_args: list, input_type: str = "contigs"
):
    """Runs the indexing using bowtie2

    Constructs and runs the final bowtie2-build command.

    Args:
        fasta_fps (list): List of FASTA files to be indexed.
        result_fp (str): Path to the result file where the indices
            will be created.
        common_args (list): List of common flags and their values for
            the bowtie2-build command.
        input_type (str): Type of input sequences. Can be mags or contigs.
    """
    _assert_inputs_not_empty(fasta_fps)

    base_cmd = ["bowtie2-build"]
    base_cmd.extend(common_args)

    for fp in fasta_fps:

        sample_dp = os.path.join(result_fp, _get_subdir_from_path(fp, input_type))
        os.makedirs(sample_dp)

        cmd = deepcopy(base_cmd)
        cmd.extend([fp, os.path.join(sample_dp, "index")])

        try:
            run_command(cmd, verbose=True)
        except subprocess.CalledProcessError as e:
            raise Exception(
                "An error was encountered while running Bowtie2, "
                f"(return code {e.returncode}), please inspect "
                "stdout and stderr to learn more."
            )


def index_contigs(
    contigs: ContigSequencesDirFmt,
    large_index: bool = False,
    debug: bool = False,
    sanitized: bool = False,
    verbose: bool = False,
    noauto: bool = False,
    packed: bool = False,
    bmax: int = "auto",
    bmaxdivn: int = 4,
    dcv: int = 1024,
    nodc: bool = False,
    offrate: int = 5,
    ftabchars: int = 10,
    threads: int = 1,
    seed: int = 0,
) -> Bowtie2IndexDirFmt:
    if bmax == "auto":
        bmax = None
    kwargs = {k: v for k, v in locals().items() if k not in ["contigs"]}
    common_args = _process_common_input_params(
        processing_func=_process_bowtie2build_arg, params=kwargs
    )
    result = Bowtie2IndexDirFmt()

    contig_fps = sorted(glob.glob(os.path.join(str(contigs), "*_contigs.fa")))
    _index_seqs(contig_fps, str(result), common_args, "contigs")

    return result


def index_mags(
    mags: MultiMAGSequencesDirFmt,
    large_index: bool = False,
    debug: bool = False,
    sanitized: bool = False,
    verbose: bool = False,
    noauto: bool = False,
    packed: bool = False,
    bmax: int = "auto",
    bmaxdivn: int = 4,
    dcv: int = 1024,
    nodc: bool = False,
    offrate: int = 5,
    ftabchars: int = 10,
    threads: int = 1,
    seed: int = 0,
) -> MultiBowtie2IndexDirFmt:
    if bmax == "auto":
        bmax = None
    kwargs = {k: v for k, v in locals().items() if k not in ["mags"]}
    common_args = _process_common_input_params(
        processing_func=_process_bowtie2build_arg, params=kwargs
    )

    result = MultiBowtie2IndexDirFmt()

    mag_fps = sorted(glob.glob(os.path.join(str(mags), "*", "*.fa*")))
    _index_seqs(mag_fps, str(result), common_args, "mags")

    return result
