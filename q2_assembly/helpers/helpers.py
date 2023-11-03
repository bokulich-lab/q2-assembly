# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import warnings

import numpy as np
from q2_types.bowtie2 import Bowtie2IndexDirFmt
from q2_types_genomics.per_sample_data import BAMDirFmt, ContigSequencesDirFmt
from qiime2.util import duplicate


def partition_contigs(
    contigs: ContigSequencesDirFmt, num_partitions: int = None
) -> ContigSequencesDirFmt:
    partitioned_contigs = {}

    contigs = list(contigs.path.iterdir())
    num_samples = len(contigs)
    if num_partitions is None:
        num_partitions = num_samples
    elif num_partitions > num_samples:
        warnings.warn(
            "You have requested a number of partitions"
            f" '{num_partitions}' that is greater than your number"
            f" of samples '{num_samples}.' Your data will be"
            f" partitioned by sample into '{num_samples}'"
            " partitions."
        )

    # If the number of splits to make is larger than the collection given,
    # np.array_split partitions the collection into individual elements
    contigs = np.array_split(contigs, num_partitions)
    for i, samples in enumerate(contigs, 1):
        result = ContigSequencesDirFmt()

        for sample_fp in samples:
            # These paths are defined in the ContigSequencesDirFmt class as
            # {sample_id}_contigs.(fa | fasta). This should get the id from a
            # name like that
            sample_id = sample_fp.name.rsplit("_contigs", 1)[0]
            duplicate(sample_fp, result.path / sample_fp.name)

        if num_partitions == num_samples:
            partitioned_contigs[sample_id] = result
        else:
            partitioned_contigs[i] = result

    return partitioned_contigs


def collate_contigs(contigs: ContigSequencesDirFmt) -> ContigSequencesDirFmt:
    collated_contigs = ContigSequencesDirFmt()

    for contig in contigs:
        for fp in contig.path.iterdir():
            duplicate(fp, collated_contigs.path / fp.name)

    return collated_contigs


def collate_indices(indices: Bowtie2IndexDirFmt) -> Bowtie2IndexDirFmt:
    collated_indices = Bowtie2IndexDirFmt()

    for index in indices:
        for in_sample_dir in index.path.iterdir():
            out_sample_dir = collated_indices.path / in_sample_dir.name
            os.mkdir(out_sample_dir)

            for fp in in_sample_dir.iterdir():
                out_fp = os.path.join(out_sample_dir, fp.name)
                duplicate(fp, out_fp)

    return collated_indices


def collate_alignments(alignments: BAMDirFmt) -> BAMDirFmt:
    collated_alignments = BAMDirFmt()

    for alignment in alignments:
        for _alignment in alignment.path.iterdir():
            filename = os.path.basename(_alignment)
            duplicate(_alignment, os.path.join(collated_alignments.path, filename))

    return collated_alignments


def _partition_sample_data(sample_data, num_partitions):
    """This one is not actually a QIIME 2 action, but it is in the vein of the
    actions here.
    """
    result = []

    num_samples = len(sample_data)
    if num_partitions is None:
        num_partitions = num_samples
    elif num_partitions > num_samples:
        warnings.warn(
            "You have requested a number of partitions"
            f" '{num_partitions}' that is greater than your number"
            f" of samples '{num_samples}.' Your data will be"
            f" partitioned by sample into '{num_samples}'"
            " partitions."
        )

    partitioned = np.array_split(list(sample_data.items()), num_partitions)
    for partition in partitioned:
        names = []
        forwards = []
        reverses = []
        indices = []
        for name, info in partition:
            names.append(name)
            forwards.append(info["fwd"])
            reverses.append(info["rev"])
            indices.append(info["index"])
        result.append((names, forwards, reverses, indices))

    return result