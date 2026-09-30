# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import contextlib
import os
import re
import shutil
import tempfile
import unittest
import uuid
from pathlib import Path
from typing import Counter
from unittest.mock import ANY, call, patch

import shortuuid
import skbio
from parameterized import parameterized
from q2_types.per_sample_sequences import AlignmentMap, BAMDirFmt, ContigSequencesDirFmt
from q2_types.sample_data import SampleData
from qiime2 import Artifact
from qiime2.plugin import Properties
from qiime2.plugin.testing import TestPluginBase

from q2_assembly.helpers.helpers import rename_contigs, sort_alignment_maps


class TestUtils(TestPluginBase):
    package = "q2_assembly.tests"
    UUID4_REGEX = re.compile(
        r"^[a-f0-9]{8}-[a-f0-9]{4}-4[a-f0-9]{3}-[89ab][a-f0-9]{3}-[a-f0-9]{12}$",
        re.IGNORECASE,
    )
    UUID3_5_REGEX = re.compile(
        r"^[a-f0-9]{8}-[a-f0-9]{4}-(3|5)[a-f0-9]{3}-[89ab][a-f0-9]{3}-[a-f0-9]{12}$",
        re.IGNORECASE,
    )

    def setUp(self):
        super().setUp()
        with contextlib.ExitStack() as stack:
            self._tmp = stack.enter_context(tempfile.TemporaryDirectory())
            self.addCleanup(stack.pop_all().close)

        self.namespace = uuid.NAMESPACE_OID
        self.name = "test-test"

    def is_valid_shortuuid(self, shortuuid_str):
        try:
            shortuuid.decode(shortuuid_str)
            return True
        except ValueError:
            return False

    def test_is_valid_shortuuid(self):
        true_shortuuid = shortuuid.uuid()
        false_shortuuid = "1234567890"
        self.assertTrue(self.is_valid_shortuuid(true_shortuuid))
        self.assertFalse(self.is_valid_shortuuid(false_shortuuid))

    @parameterized.expand(
        [
            ("uuid4", UUID4_REGEX),
            ("uuid3", UUID3_5_REGEX),
            ("uuid5", UUID3_5_REGEX),
            ("shortuuid", None),
        ]
    )
    def test_rename_contigs(self, uuid_type, regex):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")

        with tempfile.TemporaryDirectory() as tmp:
            for contig_fp in contigs.sample_dict().values():
                shutil.copyfile(
                    contig_fp, os.path.join(tmp, os.path.basename(contig_fp))
                )

            contigs_test = ContigSequencesDirFmt(tmp, "r")

            renamed_contigs = rename_contigs(
                contigs_test, uuid_type, include_sample_id=False
            )

            new_contig_ids = {
                record.metadata["id"]
                for sample_fp in renamed_contigs.sample_dict().values()
                for record in skbio.read(sample_fp, format="fasta")
            }

            # ensure the IDs are unique across samples
            # there are 14 contigs in the test data
            self.assertEqual(len(new_contig_ids), 14)

            # check if type of generated id is correct
            if uuid_type == "shortuuid":
                self.assertTrue(
                    all(self.is_valid_shortuuid(new_id) for new_id in new_contig_ids)
                )
            else:
                self.assertTrue(all(regex.match(new_id) for new_id in new_contig_ids))

    def test_rename_contigs_with_separator(self):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")

        with tempfile.TemporaryDirectory() as tmp:
            for contig_fp in contigs.sample_dict().values():
                shutil.copyfile(
                    contig_fp, os.path.join(tmp, os.path.basename(contig_fp))
                )

            contigs_test = ContigSequencesDirFmt(tmp, "r")

            renamed_contigs = rename_contigs(
                contigs_test, "shortuuid", include_sample_id=True, separator="-"
            )

            new_contig_ids = {
                record.metadata["id"]
                for sample_fp in renamed_contigs.sample_dict().values()
                for record in skbio.read(sample_fp, format="fasta")
            }

            obs_samples = Counter(name.split("-")[0] for name in new_contig_ids)
            self.assertDictEqual(obs_samples, {"sample1": 10, "sample2": 4})

    @parameterized.expand(["shortuuid", "uuid3", "uuid4", "uuid5"])
    @patch("q2_assembly.helpers.helpers.modify_contig_ids")
    def test_rename_contigs_method_call(self, uuid_type, p1):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        _ = rename_contigs(contigs, uuid_type)
        calls = []
        for sample_id, contig_fp in contigs.sample_dict().items():
            calls.append(call(ANY, sample_id, uuid_type, ":"))

        p1.assert_has_calls(calls)

    @patch("q2_assembly.helpers.helpers.run_command")
    def test_sort_alignment_maps(self, p1):
        maps = self.get_data_path("alignment_map", "r")

        with tempfile.TemporaryDirectory() as tmp:
            for map_fp in Path(maps).glob("*.bam"):
                shutil.copyfile(map_fp, os.path.join(tmp, os.path.basename(map_fp)))

            bam_dir = BAMDirFmt(tmp, "r")
            out_dir = sort_alignment_maps(bam_dir)

            calls = []
            for map_fp in Path(tmp).glob("*.bam"):
                samp_name = map_fp.stem
                sorted_bam = os.path.join(str(out_dir), f"{samp_name}.bam")
                calls.append(
                    call(
                        ["samtools", "sort", str(map_fp), "-o", sorted_bam],
                        verbose=True,
                    )
                )

            p1.assert_has_calls(calls, any_order=True)
            self.assertTrue(os.path.exists(str(out_dir)))

    @parameterized.expand(
        [
            ("dereplicated", ("mags", "dereplicated")),
            ("pooled", ("contigs", "pooled")),
            ("combined", ("contigs", "mags", "coassembled", "dereplicated")),
        ]
    )
    @patch("q2_assembly.helpers.helpers.run_command")
    def test_sort_and_collate_preserve_properties(self, name, properties, command):
        # Materialize a valid output while keeping samtools outside this unit test.
        command.side_effect = lambda args, **kwargs: shutil.copyfile(args[2], args[4])
        maps = Artifact.import_data(
            SampleData[AlignmentMap % Properties(properties)],
            self.get_data_path("alignment_map"),
        )
        collate = self.plugin.methods["collate_alignments"]
        sort = self.plugin.methods["sort_alignment_maps"]

        (unsorted,) = collate([maps])
        self.assertEqual(unsorted.type, maps.type)
        (sorted_maps,) = sort(unsorted)
        expected = SampleData[AlignmentMap % Properties((*properties, "sorted"))]
        self.assertEqual(sorted_maps.type, expected)
        (collated,) = collate([sorted_maps])
        self.assertEqual(collated.type, expected)
        collated.validate()


if __name__ == "__main__":
    unittest.main()
