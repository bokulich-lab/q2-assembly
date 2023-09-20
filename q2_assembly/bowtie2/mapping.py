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
from copy import deepcopy

import pandas as pd
from q2_types.bowtie2 import Bowtie2IndexDirFmt
from q2_types.per_sample_sequences import SingleLanePerSamplePairedEndFastqDirFmt
from q2_types_genomics.per_sample_data import BAMDirFmt

from q2_assembly.helpers.helpers import _partition_sample_data

from .._utils import _process_common_input_params, run_commands_with_pipe
from .utils import _process_bowtie2_arg


def map_reads_to_contigs(
    ctx,
    indexed_contigs,
    reads,
    skip=0,
    qupto="unlimited",
    trim5=0,
    trim3=0,
    trim_to="untrimmed",
    phred33=False,
    phred64=False,
    mode="local",
    sensitivity="sensitive",
    n=0,
    len=22,
    i="S,1,1.15",
    n_ceil="L,0,0.15",
    dpad=15,
    gbar=4,
    ignore_quals=False,
    nofw=False,
    norc=False,
    no_1mm_upfront=False,
    end_to_end=False,
    local=False,
    ma=2,
    mp=6,
    np=1,
    rdg="5,3",
    rfg="5,3",
    k="off",
    a=False,
    d=15,
    r=2,
    minins=0,
    maxins=500,
    valid_mate_orientations="fr",
    no_mixed=False,
    no_discordant=False,
    dovetail=False,
    no_contain=False,
    no_overlap=False,
    offrate="off",
    threads=1,
    reorder=False,
    mm=False,
    seed=0,
    non_deterministic=False,
    num_partitions=None,
):
    if qupto == "unlimited":
        qupto = None
    if trim_to == "untrimmed":
        trim_to = None
    if k == "off":
        k = None
    if offrate == "off":
        offrate = None
    kwargs = {
        k: v
        for k, v in locals().items()
        if k
        not in [
            "ctx",
            "indexed_contigs",
            "reads",
            "sensitivity",
            "mode",
            "num_partitions",
        ]
    }

    _map_sample_reads = ctx.get_action("assembly", "_map_sample_reads")
    collate_alignments = ctx.get_action("assembly", "collate_alignments")

    common_args = _process_common_input_params(
        processing_func=_process_bowtie2_arg, params=kwargs
    )
    if mode == "local":
        common_args.append(f"--{sensitivity}-{mode}")
    else:
        common_args.append(f"--{sensitivity}")

    indexed_contigs_format = indexed_contigs.view(indexed_contigs.format)
    reads_manifest = reads.view(reads.format).manifest.view(pd.DataFrame)
    paired = isinstance(reads.format, SingleLanePerSamplePairedEndFastqDirFmt)
    sample_data = _gather_sample_data(indexed_contigs_format, reads_manifest, paired)

    mapped_reads = []
    partitioned_sample_data = _partition_sample_data(sample_data, num_partitions)
    for partition in partitioned_sample_data:
        names, forwards, reverses, indices = partition
        (mapped_read,) = _map_sample_reads(
            common_args, paired, names, forwards, reverses, indices
        )
        mapped_reads.append(mapped_read)

    (collated_mapped_reads,) = collate_alignments(mapped_reads)
    return collated_mapped_reads


def _gather_sample_data(
    indexed_contigs: Bowtie2IndexDirFmt, reads_manifest: pd.DataFrame, paired: bool
):
    """Collects reads and indices for all the samples.

    Args:
        indexed_contigs (Bowtie2IndexDirFmt): Indices for all contigs.
        reads_manifest (pd.DataFrame): Manifest with forward and
            reverse (if paired) reads.
        paired (bool): Indicates whether reads are paired-end.

    Returns:
        sample_data (dict): Dictionary with read and index information per sample.
    """
    full_set = {}
    for samp in list(reads_manifest.index):
        full_set[samp] = {
            "fwd": reads_manifest.loc[samp, "forward"],
            "rev": reads_manifest.loc[samp, "reverse"] if paired else "",
        }
    for samp in full_set.keys():
        indpath = os.path.join(str(indexed_contigs), samp)
        if os.path.exists(indpath):
            full_set[samp]["index"] = os.path.join(indpath, "index")
        else:
            raise Exception(
                f"Index files missing for sample {samp}. " f"Please check your input."
            )
    return full_set


def _map_sample_reads(
    common_args: list,
    paired: bool,
    sample_names: list,
    sample_fwds: list,
    sample_revs: list,
    sample_indices: list,
) -> BAMDirFmt:
    """Constructs mapping command for bowtie2 and runs the mapping.

    Args:
        common_args (list): List of arguments that should be passed to bowtie2.
        paired (bool): Indicates whether reads are paired-end.
        sample_name (str): Name of the sample to be processed.
        sample_inputs (dict): Dictionary of sample inputs containing
            bowtie2 index, forward and reverse (if paired) reads.
        result_fp (str): Path to where the result should be stored.
    """
    result = BAMDirFmt()

    for name, fwd, rev, index in zip(
        sample_names, sample_fwds, sample_revs, sample_indices
    ):
        base_cmd = ["bowtie2"]
        base_cmd.extend(common_args)

        with tempfile.TemporaryDirectory() as tmp:
            bam_result = os.path.join(tmp, "alignment.bam")
            cmd = deepcopy(base_cmd)
            cmd.extend(
                [
                    "-x",
                    index,
                ]
            )

            if paired:
                cmd.extend(["-1", fwd, "-2", rev])
            else:
                cmd.extend(["-U", fwd])

            try:
                run_commands_with_pipe(
                    cmd1=cmd,
                    cmd2=["samtools", "view", "-bS", "-o", bam_result],
                    verbose=True,
                )
            except subprocess.CalledProcessError as e:
                raise Exception(
                    "An error was encountered while running Bowtie2, "
                    f"(return code {e.returncode}), please inspect "
                    "stdout and stderr to learn more."
                )

            shutil.move(bam_result, os.path.join(str(result), f"{name}_alignment.bam"))

            # TODO: introduce support for the stats generated by bowtie
            #  (--met flag) and create a visualisation

    return result
