# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
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
import warnings
from unittest.mock import ANY, call, patch

import shortuuid
import skbio
from parameterized import parameterized
from q2_types.feature_data import DNAFASTAFormat
from q2_types.genome_data import GenomeSequencesDirectoryFormat
from q2_types.per_sample_sequences import ContigSequencesDirFmt
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import assembly

from q2_assembly.helpers.helpers import collate_genomes, rename_contigs


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

    @parameterized.expand(["single", "multiple"])
    def test_collate_genomes_dnafastaformat(self, input):
        genomes1 = DNAFASTAFormat(
            self.get_data_path("dna-fasta-format/dna-sequences1.fasta"), "r"
        )
        genomes2 = DNAFASTAFormat(
            self.get_data_path("dna-fasta-format/dna-sequences2.fasta"), "r"
        )
        if input == "single":
            genomes = [genomes1]
            content = {
                "ref1": {"description": "d_Bacteria_1", "sequence": "ACGTACGT"},
                "ref2": {"description": "d_Bacteria_2", "sequence": "CGTCGTCC"},
            }
            exp_files = ["ref1.fasta", "ref2.fasta"]
        else:
            genomes = [genomes1, genomes2]
            content = {
                "ref1": {"description": "d_Bacteria_1", "sequence": "ACGTACGT"},
                "ref2": {"description": "d_Bacteria_2", "sequence": "CGTCGTCC"},
                "ref5": {"description": "d_Bacteria_3", "sequence": "ACGTACGT"},
                "ref6": {"description": "d_Bacteria_4", "sequence": "CGTCGTCC"},
            }
            exp_files = ["ref1.fasta", "ref2.fasta", "ref5.fasta", "ref6.fasta"]

        collated_genomes = collate_genomes(genomes=genomes)
        actual_files = sorted(os.listdir(collated_genomes.path))
        self.assertEqual(actual_files, exp_files)

        for fn in actual_files:
            fp = os.path.join(collated_genomes.path, fn)
            with open(fp, "r") as fasta_file:
                for seq in skbio.io.read(fasta_file, "fasta"):
                    actual_id = seq.metadata["id"]
                    actual_description = seq.metadata["description"]
                    actual_sequence = str(seq)
                    expected_id = fn.split(".")[0]
                    expected_desc = content[expected_id]["description"]
                    expected_sequence = content[expected_id]["sequence"]

                    self.assertEquals(actual_id, expected_id)
                    self.assertEqual(actual_description, expected_desc)
                    self.assertEqual(actual_sequence, expected_sequence)

    def test_collate_genomes_genome_dir_multiple(self):
        genomes1 = GenomeSequencesDirectoryFormat(
            self.get_data_path("genomes-dir-format1"), "r"
        )
        genomes2 = GenomeSequencesDirectoryFormat(
            self.get_data_path("genomes-dir-format2"), "r"
        )
        genomes = [genomes1, genomes2]
        collated_genomes = collate_genomes(genomes=genomes)
        exp_files = ["ref1.fasta", "ref2.fasta", "ref3.fasta"]
        actual_files = sorted(os.listdir(collated_genomes.path))
        self.assertEqual(exp_files, actual_files)

    def test_collate_genomes_mix(self):
        # should throw TypeError
        genomes1 = DNAFASTAFormat(
            self.get_data_path("dna-fasta-format/dna-sequences1.fasta"), "r"
        )
        genomes2 = GenomeSequencesDirectoryFormat(
            self.get_data_path("genomes-dir-format2"), "r"
        )
        genomes = [genomes2, genomes1]
        with self.assertRaises(TypeError):
            assembly.methods.collate_genomes(genomes=genomes)

    @parameterized.expand(["GenomeData", "DNAFASTAFormat"])
    def test_collate_genomes_dnafastaformat_multiple_duplicates_warn(self, dir_fmt):
        duplicate_ids = (
            ["ref1.fasta", "ref2.fasta"]
            if dir_fmt == "GenomeData"
            else ["ref1", "ref2"]
        )
        warn_msg = (
            "Duplicate sequence files were found for the following IDs: {}. "
            "The latest occurrence will overwrite all previous occurrences "
            "for each corresponding ID."
        ).format(", ".join(duplicate_ids))
        if dir_fmt == "GenomeData":
            genomes1 = GenomeSequencesDirectoryFormat(
                self.get_data_path("genomes-dir-format1"), "r"
            )
        else:
            genomes1 = DNAFASTAFormat(
                self.get_data_path("dna-fasta-format/dna-sequences1.fasta"), "r"
            )
        with warnings.catch_warnings(record=True) as w:
            collated_genomes = collate_genomes(genomes=[genomes1, genomes1])
            exp_files = ["ref1.fasta", "ref2.fasta"]
            actual_files = sorted(os.listdir(collated_genomes.path))
            self.assertEqual(actual_files, exp_files)
            self.assertEqual(warn_msg, str(w[0].message))

            if dir_fmt == "DNAFASTAFormat":
                content = {
                    "ref1": {"description": "d_Bacteria_1", "sequence": "ACGTACGT"},
                    "ref2": {"description": "d_Bacteria_2", "sequence": "CGTCGTCC"},
                }

                for fn in actual_files:
                    fp = os.path.join(collated_genomes.path, fn)
                    with open(fp, "r") as fasta_file:
                        for seq in skbio.io.read(fasta_file, "fasta"):
                            actual_id = seq.metadata["id"]
                            actual_description = seq.metadata["description"]
                            actual_sequence = str(seq)
                            expected_id = fn.split(".")[0]
                            expected_desc = content[expected_id]["description"]
                            expected_sequence = content[expected_id]["sequence"]

                            self.assertEquals(actual_id, expected_id)
                            self.assertEqual(actual_description, expected_desc)
                            self.assertEqual(actual_sequence, expected_sequence)

    @parameterized.expand(["GenomeData", "DNAFASTAFormat"])
    def test_collate_genomes_duplicates_error(self, dir_fmt):
        duplicate_ids = ["ref3.fasta"] if dir_fmt == "GenomeData" else ["ref1"]
        error_msg = (
            "Duplicate sequence files were found for the "
            "following IDs: %s." % ", ".join(duplicate_ids)
        )
        if dir_fmt == "GenomeData":
            genomes1 = GenomeSequencesDirectoryFormat(
                self.get_data_path("genomes-dir-format2"), "r"
            )
        else:
            genomes1 = DNAFASTAFormat(
                self.get_data_path("dna-fasta-format/dna-sequences1.fasta"), "r"
            )
        with self.assertRaisesRegex(ValueError, error_msg):
            _ = collate_genomes(genomes=[genomes1, genomes1], on_duplicates="error")

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

            renamed_contigs = rename_contigs(contigs_test, uuid_type)

            new_contig_ids = set()

            i = 0
            for sample_id, sample_fp in renamed_contigs.sample_dict().items():
                for record in skbio.read(sample_fp, format="fasta"):
                    new_contig_ids.add(record.metadata["id"])
                    i = i + 1

            self.assertEqual(
                len(new_contig_ids), i
            )  # ensure the IDs are unique across samples

            # check if type of generated id is correct
            if uuid_type == "shortuuid":
                self.assertTrue(
                    all(self.is_valid_shortuuid(new_id) for new_id in new_contig_ids)
                )
            else:
                self.assertTrue(all(regex.match(new_id) for new_id in new_contig_ids))

    @parameterized.expand(["shortuuid", "uuid3", "uuid4", "uuid5"])
    @patch("q2_assembly.helpers.modify_contig_ids")
    def test_rename_contigs_method_call(self, uuid_type, p1):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        _ = rename_contigs(contigs, uuid_type)
        calls = []
        for sample_id, contig_fp in contigs.sample_dict().items():
            calls.append(call(ANY, sample_id, uuid_type))

        p1.assert_has_calls(calls)


if __name__ == "__main__":
    unittest.main()
