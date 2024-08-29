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
from unittest.mock import ANY, call, patch

import shortuuid
import skbio
from parameterized import parameterized
from q2_types._util import DNAFASTAFormat
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

    def test_collate_genomes_dnafastaformat_single(self):
        genomes1 = DNAFASTAFormat(
            self.get_data_path("dna-fasta-format/dna-sequences1.fasta"), "r"
        )
        collated_genomes = collate_genomes(genomes_in=[genomes1])
        self.assertEqual(len(os.listdir(collated_genomes.path)), 2)

    def test_collate_genomes_dnafastaformat_multiple(self):
        genomes1 = DNAFASTAFormat(
            self.get_data_path("dna-fasta-format/dna-sequences1.fasta"), "r"
        )
        genomes2 = DNAFASTAFormat(
            self.get_data_path("dna-fasta-format/dna-sequences2.fasta"), "r"
        )
        collated_genomes = collate_genomes(genomes_in=[genomes1, genomes2])
        self.assertEqual(len(os.listdir(collated_genomes.path)), 4)

    def test_collate_genomes_genome_dir_multiple(self):
        genomes1 = GenomeSequencesDirectoryFormat(
            self.get_data_path("genomes-dir-format1"), "r"
        )
        genomes2 = GenomeSequencesDirectoryFormat(
            self.get_data_path("genomes-dir-format2"), "r"
        )
        genomes = [genomes1, genomes2]
        collated_genomes = collate_genomes(genomes_in=genomes)
        self.assertEqual(len(os.listdir(collated_genomes.path)), 4)

    def test_collate_genomes_mix(self):
        # should throw TypeError
        genomes1 = DNAFASTAFormat(
            self.get_data_path("dna-fasta-format/dna-sequences1.fasta"), "r"
        )
        genomes2 = GenomeSequencesDirectoryFormat(
            self.get_data_path("genomes-dir-format1"), "r"
        )

        genomes = [genomes1, genomes2]

        with self.assertRaises(TypeError):
            # collate_genomes(genomes_in=genomes) this does not throw the exception
            assembly.methods.collate_genomes(genomes=genomes)

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
