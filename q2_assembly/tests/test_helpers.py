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

from q2_assembly._utils import generate_shortuuid
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

    @patch("q2_assembly._utils._generate_unique_uuid")
    def test_rename_contigs(self, p1):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        p1.return_value = generate_shortuuid()
        renamed_contigs = rename_contigs(contigs, "shortuuid")
        set_of_contig_ids = set()

        for sample_id, sample_fp in renamed_contigs.sample_dict().items():
            for record in skbio.read(sample_fp, format="fasta"):
                set_of_contig_ids.add(record.metadata["id"])

        self.assertEqual(len(set_of_contig_ids), 1)
        self.assertEqual(list(set_of_contig_ids)[0], p1.return_value)

    @patch("q2_assembly.helpers.modify_contig_ids")
    def test_rename_contigs_method_call(self, p1):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        _ = rename_contigs(contigs, "uuid4")
        p1.assert_has_calls([call(ANY, ANY, "uuid4"), call(ANY, ANY, "uuid4")])


if __name__ == "__main__":
    unittest.main()
