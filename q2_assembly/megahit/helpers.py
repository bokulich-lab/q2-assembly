# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types_genomics.per_sample_data import ContigSequencesDirFmt
from qiime2.util import duplicate


def partition_contigs(
    contigs: ContigSequencesDirFmt, num_partitions: int = None
) -> ContigSequencesDirFmt:
    collection = {}
    for sample_fp in contigs.path.iterdir():
        result = ContigSequencesDirFmt()
        # These paths are defined in the ContigSequencesDirFmt class as
        # {sample_id}_contigs,(fa | fasta). This should get the id from a name
        # like that
        sample_id = sample_fp.name.split("_contigs")[0]
        duplicate(sample_fp, result.path / sample_fp.name)

        collection[sample_id] = result

    return collection


def collate_contigs(contigs: ContigSequencesDirFmt) -> ContigSequencesDirFmt:
    collated_contigs = ContigSequencesDirFmt()

    for contig in contigs:
        for fp in contig.path.iterdir():
            duplicate(fp, collated_contigs.path / fp.name)

    return collated_contigs
