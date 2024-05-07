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
from typing import Union

import pandas as pd
from q2_types.bowtie2 import Bowtie2IndexDirFmt
from q2_types.feature_data import FeatureData
from q2_types.per_sample_sequences import (
    BAMDirFmt,
    MultiBowtie2Index,
    MultiBowtie2IndexDirFmt,
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
    SingleBowtie2Index,
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from q2_types.sample_data import SampleData
from qiime2.core.type import Properties

from .._utils import _process_common_input_params, run_commands_with_pipe
from .utils import _process_bowtie2_arg


def map_reads(
    ctx,
    index,
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
    kwargs = {
        k: v
        for k, v in locals().items()
        if k
        not in [
            "ctx",
            "reads",
            "num_partitions",
        ]
    }

    collate_alignments = ctx.get_action("assembly", "collate_alignments")

    if reads.type <= SampleData[SequencesWithQuality]:
        partition_method = ctx.get_action("demux", "partition_samples_single")
    elif reads.type <= SampleData[PairedEndSequencesWithQuality]:
        partition_method = ctx.get_action("demux", "partition_samples_paired")
    else:
        raise NotImplementedError()

    if index.type <= SampleData[SingleBowtie2Index % Properties("contigs")]:
        _map_reads = ctx.get_action("assembly", "_map_reads_to_contigs")
    elif index.type <= SampleData[MultiBowtie2Index % Properties("mags")]:
        _map_reads = ctx.get_action("assembly", "_map_reads_to_mags")
    elif index.type <= FeatureData[SingleBowtie2Index % Properties("mags")]:
        _map_reads = ctx.get_action("assembly", "_map_reads_to_mags")
    else:
        raise NotImplementedError()

    (partitioned_reads,) = partition_method(reads, num_partitions)

    mapped_reads = []
    for read in partitioned_reads.values():
        (mapped_read,) = _map_reads(reads=read, **kwargs)
        mapped_reads.append(mapped_read)

    (collated_mapped_reads,) = collate_alignments(mapped_reads)
    return collated_mapped_reads


def _map_reads_to_contigs(
    index: Bowtie2IndexDirFmt,
    reads: Union[
        SingleLanePerSamplePairedEndFastqDirFmt, SingleLanePerSampleSingleEndFastqDirFmt
    ],
    skip: int = 0,
    qupto: int = "unlimited",
    trim5: int = 0,
    trim3: int = 0,
    trim_to: str = "untrimmed",
    phred33: bool = False,
    phred64: bool = False,
    mode: str = "local",
    sensitivity: str = "sensitive",
    n: int = 0,
    len: int = 22,
    i: str = "S,1,1.15",
    n_ceil: str = "L,0,0.15",
    dpad: int = 15,
    gbar: int = 4,
    ignore_quals: bool = False,
    nofw: bool = False,
    norc: bool = False,
    no_1mm_upfront: bool = False,
    end_to_end: bool = False,
    local: bool = False,
    ma: int = 2,
    mp: int = 6,
    np: int = 1,
    rdg: str = "5,3",
    rfg: str = "5,3",
    k: int = "off",
    a: bool = False,
    d: int = 15,
    r: int = 2,
    minins: int = 0,
    maxins: int = 500,
    valid_mate_orientations: str = "fr",
    no_mixed: bool = False,
    no_discordant: bool = False,
    dovetail: bool = False,
    no_contain: bool = False,
    no_overlap: bool = False,
    offrate: int = "off",
    threads: int = 1,
    reorder: bool = False,
    mm: bool = False,
    seed: int = 0,
    non_deterministic: bool = False,
) -> BAMDirFmt:
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
        if k not in ["index", "reads", "sensitivity", "mode"]
    }
    common_args = _process_common_input_params(
        processing_func=_process_bowtie2_arg, params=kwargs
    )
    if mode == "local":
        common_args.append(f"--{sensitivity}-{mode}")
    else:
        common_args.append(f"--{sensitivity}")

    paired = isinstance(reads, SingleLanePerSamplePairedEndFastqDirFmt)
    manifest = reads.manifest.view(pd.DataFrame)

    full_sample_set = _gather_sample_data(index, manifest, paired)

    result = BAMDirFmt()
    for samp, props in full_sample_set.items():
        _map_sample_reads(common_args, paired, samp, props, str(result))

    return result


def _map_reads_to_mags(
    index: Bowtie2IndexDirFmt,
    reads: Union[
        SingleLanePerSamplePairedEndFastqDirFmt, SingleLanePerSampleSingleEndFastqDirFmt
    ],
    skip: int = 0,
    qupto: int = "unlimited",
    trim5: int = 0,
    trim3: int = 0,
    trim_to: str = "untrimmed",
    phred33: bool = False,
    phred64: bool = False,
    mode: str = "local",
    sensitivity: str = "sensitive",
    n: int = 0,
    len: int = 22,
    i: str = "S,1,1.15",
    n_ceil: str = "L,0,0.15",
    dpad: int = 15,
    gbar: int = 4,
    ignore_quals: bool = False,
    nofw: bool = False,
    norc: bool = False,
    no_1mm_upfront: bool = False,
    end_to_end: bool = False,
    local: bool = False,
    ma: int = 2,
    mp: int = 6,
    np: int = 1,
    rdg: str = "5,3",
    rfg: str = "5,3",
    k: int = "off",
    a: bool = False,
    d: int = 15,
    r: int = 2,
    minins: int = 0,
    maxins: int = 500,
    valid_mate_orientations: str = "fr",
    no_mixed: bool = False,
    no_discordant: bool = False,
    dovetail: bool = False,
    no_contain: bool = False,
    no_overlap: bool = False,
    offrate: int = "off",
    threads: int = 1,
    reorder: bool = False,
    mm: bool = False,
    seed: int = 0,
    non_deterministic: bool = False,
) -> BAMDirFmt:
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
        if k not in ["index", "reads", "sensitivity", "mode"]
    }
    common_args = _process_common_input_params(
        processing_func=_process_bowtie2_arg, params=kwargs
    )
    if mode == "local":
        common_args.append(f"--{sensitivity}-{mode}")
    else:
        common_args.append(f"--{sensitivity}")

    paired = isinstance(reads, SingleLanePerSamplePairedEndFastqDirFmt)
    manifest = reads.manifest.view(pd.DataFrame)

    full_sample_set = _gather_feature_data(index, manifest, paired)

    result = BAMDirFmt()
    for samp, props in full_sample_set.items():
        _map_sample_reads(common_args, paired, samp, props, str(result))

    return result


def _map_sample_reads(
    common_args: list,
    paired: bool,
    sample_name: str,
    sample_inputs: dict,
    result_fp: str,
):
    """Constructs mapping command for bowtie2 and runs the mapping.

    Args:
        common_args (list): List of arguments that should be passed to bowtie2.
        paired (bool): Indicates whether reads are paired-end.
        sample_name (str): Name of the sample to be processed.
        sample_inputs (dict): Dictionary of sample inputs containing
            bowtie2 index, forward and reverse (if paired) reads.
        result_fp (str): Path to where the result should be stored.
    """
    base_cmd = ["bowtie2"]
    base_cmd.extend(common_args)

    with tempfile.TemporaryDirectory() as tmp:
        bam_result = os.path.join(tmp, "alignment.bam")
        cmd = deepcopy(base_cmd)
        cmd.extend(
            [
                "-x",
                sample_inputs["index"],
            ]
        )

        if paired:
            cmd.extend(["-1", sample_inputs["fwd"], "-2", sample_inputs["rev"]])
        else:
            cmd.extend(["-U", sample_inputs["fwd"]])

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

        shutil.move(bam_result, os.path.join(result_fp, f"{sample_name}_alignment.bam"))

        # TODO: introduce support for the stats generated by bowtie
        #  (--met flag) and create a visualisation


def _gather_sample_data(
    index: Union[Bowtie2IndexDirFmt, MultiBowtie2IndexDirFmt],
    reads_manifest: pd.DataFrame,
    paired: bool,
):
    """Collects reads and indices for all the samples.

    Bowtie2IndexDirFmt represents contig indices per sample (generated
    from SampleData[Contigs]) or MAG indices per sample (generated
    from SampleData[MAGs]), while MultiBowtie2IndexDirFmt represents
    MAG indices per sample per MAG (also generated from SampleData[MAGs]).
    In all cases mapping would be performed per sample.

    Args:
        index (Bowtie2IndexDirFmt | MultiBowtie2IndexDirFmt):
            Indices for all contigs/MAGs.
        reads_manifest (pd.DataFrame): Manifest with forward and
            reverse (if paired) reads.
        paired (bool): Indicates whether reads are paired-end.

    Returns:
        full_set (dict): Dictionary with read and index information per sample.
    """
    full_set = {}
    for samp in list(reads_manifest.index):
        full_set[samp] = {
            "fwd": reads_manifest.loc[samp, "forward"],
            "rev": reads_manifest.loc[samp, "reverse"] if paired else None,
        }
    for samp in full_set.keys():
        # case 1: indices per sample
        if isinstance(index, Bowtie2IndexDirFmt):
            indpath = os.path.join(str(index), samp)
        # case 2: indices per sample per MAG
        elif isinstance(index, MultiBowtie2IndexDirFmt):
            raise NotImplementedError(
                "Mapping to MAGs on a MAG-by-MAG basis is not yet supported."
            )
        else:
            raise NotImplementedError()

        if os.path.exists(indpath):
            full_set[samp]["index"] = os.path.join(indpath, "index")
        else:
            raise Exception(
                f"Index files missing for sample {samp}. " f"Please check your input."
            )
    return full_set


def _gather_feature_data(
    index: Bowtie2IndexDirFmt,
    reads_manifest: pd.DataFrame,
    paired: bool,
):
    """Collects reads and indices for all MAGs.

    We assume that all MAGs were indexed together - there will be only
    one set of Bowtie2 index files.

    Args:
        index (Bowtie2IndexDirFmt): Indices for all dereplicated MAGs.
        reads_manifest (pd.DataFrame): Manifest with forward and
            reverse (if paired) reads.
        paired (bool): Indicates whether reads are paired-end.

    Returns:
        full_set (dict): Dictionary with read and index information per sample.
    """
    full_set = {}
    for samp in list(reads_manifest.index):
        full_set[samp] = {
            "fwd": reads_manifest.loc[samp, "forward"],
            "rev": reads_manifest.loc[samp, "reverse"] if paired else None,
        }
        full_set[samp]["index"] = os.path.join(str(index), "index")

    return full_set
