# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import contextlib
import tempfile
import unittest
import uuid
from unittest.mock import ANY, call, patch

import skbio
from q2_types.per_sample_sequences import ContigSequencesDirFmt
from qiime2.plugin.testing import TestPluginBase

from q2_assembly.helpers.helpers import rename_contigs


class TestUtils(TestPluginBase):
    package = "q2_assembly.tests"

    def setUp(self):
        super().setUp()
        with contextlib.ExitStack() as stack:
            self._tmp = stack.enter_context(tempfile.TemporaryDirectory())
            self.addCleanup(stack.pop_all().close)

        self.namespace = uuid.NAMESPACE_OID
        self.name = "test-test"

    def test_rename_contigs(self):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        new_id_contigs = rename_contigs(contigs, "shortuuid")
        set_of_contig_ids = set()

        i = 0
        for sample_id, sample_fp in new_id_contigs.sample_dict().items():
            for record in skbio.read(sample_fp, format="fasta"):
                set_of_contig_ids.add(record.metadata["id"])
                i = i + 1

        self.assertEqual(
            len(set_of_contig_ids), i
        )  # ensure the IDs are unique across samples

    @patch("q2_assembly.helpers.modify_contig_ids")
    def test_rename_contigs_method_call(self, p1):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        _ = rename_contigs(contigs, "uuid4")
        p1.assert_has_calls([call(ANY, ANY, "uuid4"), call(ANY, ANY, "uuid4")])


if __name__ == "__main__":
    unittest.main()
