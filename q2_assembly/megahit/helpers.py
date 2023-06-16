# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.util import duplicate
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt


def collate_contigs(contigs: ContigSequencesDirFmt) -> ContigSequencesDirFmt:
    result = ContigSequencesDirFmt()

    for contig in contigs:
        for fp in contig.path.iterdir():
            duplicate(result.path, fp)

    return result
