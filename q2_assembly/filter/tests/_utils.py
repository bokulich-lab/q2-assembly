# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import qiime2 as q2
import skbio
from q2_types.per_sample_sequences import ContigSequencesDirFmt
from qiime2.plugin.testing import TestPluginBase


class FilterContigsTestBase(TestPluginBase):
    package = "q2_assembly.tests"

    @staticmethod
    def _read_contig_ids(fp):
        return [c.metadata["id"] for c in skbio.io.read(fp, format="fasta")]

    def _read_all_contig_ids(self, contigs):
        if isinstance(contigs, q2.Artifact):
            contigs = contigs.view(ContigSequencesDirFmt)
        return [
            _id
            for fp in contigs.sample_dict().values()
            for _id in self._read_contig_ids(fp)
        ]
