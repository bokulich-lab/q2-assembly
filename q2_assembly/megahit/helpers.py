# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types_genomics.per_sample_data import ContigSequencesDirFmt
from qiime2.util import duplicate


def collate_contigs(contigs: ContigSequencesDirFmt) -> ContigSequencesDirFmt:
    result = ContigSequencesDirFmt()

    for contig in contigs:
        for fp in contig.path.iterdir():
            duplicate(fp, result.path / fp.name)

    return result
