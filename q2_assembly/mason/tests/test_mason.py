# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import unittest
from unittest.mock import call, patch, MagicMock

from q2_types.genome_data import GenomeSequencesDirectoryFormat
from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt
from qiime2.plugin.testing import TestPluginBase

from q2_assembly.mason.mason import (
    _process_mason_arg,
    _simulate_reads_mason,
    generate_abundances,
    abundances_to_df,
    _combine_reads,
    _process_sample,
    simulate_reads_mason,
)


class MockTempDir(tempfile.TemporaryDirectory):
    pass


class TestMason(TestPluginBase):
    package = "q2_assembly.tests"

    def setUp(self):
        super().setUp()
        self.refs = GenomeSequencesDirectoryFormat(
            self.get_data_path("genomes-dir-format1"), "r"
        )

    def test_process_mason_arg(self):
        obs = _process_mason_arg("num_reads", 1000000)
        exp = ["--num_reads", "1000000"]
        self.assertListEqual(obs, exp)

    def test_process_mason_arg_list(self):
        obs = _process_mason_arg("k_list", [1, 2, 3])
        exp = ["--k_list", "1", "2", "3"]
        self.assertListEqual(obs, exp)

    def test_process_mason_arg_bool(self):
        obs = _process_mason_arg("k_bool", True)
        exp = ["--k_bool"]
        self.assertListEqual(obs, exp)

    @patch("q2_assembly.mason._process_sample")
    def test_simulate_reads_mason_helper(self, p_process):
        mock_genomes_dir_fmt = GenomeSequencesDirectoryFormat()

        refs = GenomeSequencesDirectoryFormat(
            self.get_data_path("genomes-dir-format1"), "r"
        )

        with patch(
            "q2_assembly.mason.GenomeSequencesDirectoryFormat",
            return_value=mock_genomes_dir_fmt,
        ):
            result = _simulate_reads_mason(
                reference_genomes=refs,
                sample_name="sample1",
                num_reads=2000000,
                read_length=125,
                abundance_profile="uniform",
            )
        self.assertIsInstance(result, CasavaOneEightSingleLanePerSampleDirFmt)
        p_process.assert_called_once_with(
            sample="sample1",
            genome_files=[
                os.path.join(str(mock_genomes_dir_fmt), "ref1.fasta"),
                os.path.join(str(mock_genomes_dir_fmt), "ref2.fasta"),
            ],
            abundances=[0.5, 0.5],
            total_reads=2000000,
            results_dir=str(result),
            threads=1,
            read_len=125,
            seed=42,
        )

    def test_simulate_reads_mason_one_sample(self):
        f1 = MagicMock(return_value=("reads",))
        f2 = MagicMock(return_value=("collated_reads",))
        mock_action = MagicMock(side_effect=[f1, f2])
        mock_ctx = MagicMock(get_action=mock_action)

        simulate_reads_mason(
            ctx=mock_ctx,
            reference_genomes=self.refs,
            sample_names=["sample1"],
            num_reads=[2000000],
            read_length=[125],
            abundance_profiles=["uniform"],
        )

        f1.assert_called_once_with(
            reference_genomes=self.refs,
            sample_name="sample1",
            num_reads=2000000,
            read_length=125,
            abundance_profile="uniform",
            random_seed=42,
            threads=1,
        )
        f2.assert_called_once_with(["reads"])

    def test_simulate_reads_mason_multiple_samples(self):
        f1 = MagicMock(side_effect=(("reads1",), ("reads2",)))
        f2 = MagicMock(return_value=("collated_reads",))
        mock_action = MagicMock(side_effect=[f1, f2])
        mock_ctx = MagicMock(get_action=mock_action)

        simulate_reads_mason(
            ctx=mock_ctx,
            reference_genomes=self.refs,
            sample_names=["sample1", "sample2"],
            num_reads=[2000, 4000],
            read_length=[125, 150],
            abundance_profiles=["uniform", "lognormal"],
        )

        f1.assert_has_calls(
            [
                call(
                    reference_genomes=self.refs,
                    sample_name="sample1",
                    num_reads=2000,
                    read_length=125,
                    abundance_profile="uniform",
                    random_seed=42,
                    threads=1,
                ),
                call(
                    reference_genomes=self.refs,
                    sample_name="sample2",
                    num_reads=4000,
                    read_length=150,
                    abundance_profile="lognormal",
                    random_seed=42,
                    threads=1,
                ),
            ]
        )
        f2.assert_called_once_with(["reads1", "reads2"])

    def test_simulate_reads_mason_duplicate_sample_names(self):
        with self.assertRaisesRegex(ValueError, "Sample names need to be unique"):
            simulate_reads_mason(
                ctx=MagicMock(),
                reference_genomes=self.refs,
                sample_names=["sample1", "sample1"],
                abundance_profiles=["uniform", "lognormal"],
                num_reads=[1000000, 2000000],
                read_length=[100, 150],
                threads=1,
            )

    def test_simulate_reads_mason_num_reads_length(self):
        with self.assertRaisesRegex(
            ValueError,
            "The length of read_length must be either 1 "
            "or equal to the number of sample names",
        ):
            simulate_reads_mason(
                ctx=MagicMock(),
                reference_genomes=self.refs,
                sample_names=["sample1", "sample2"],
                abundance_profiles=["uniform", "lognormal"],
                num_reads=[1000000],
                read_length=[100, 150, 200],
                threads=1,
            )

    def test_simulate_reads_mason_num_reads(self):
        with self.assertRaisesRegex(
            ValueError,
            "The length of num_reads must be either 1 "
            "or equal to the number of sample names",
        ):
            simulate_reads_mason(
                ctx=MagicMock(),
                reference_genomes=self.refs,
                sample_names=["sample1", "sample2", "sample3"],
                num_reads=[1000000, 2000000],
                read_length=[100, 150, 200],
                abundance_profiles=["uniform", "lognormal", "exponential"],
                threads=1,
            )

    def test_simulate_reads_mason_abundance_profiles(self):
        with self.assertRaisesRegex(
            ValueError,
            "The length of abundance_profiles must be either 1 "
            "or equal to the number of sample names",
        ):
            simulate_reads_mason(
                ctx=MagicMock(),
                reference_genomes=self.refs,
                sample_names=["sample1", "sample2", "sample3"],
                num_reads=[1000, 2000, 3000],
                read_length=[100, 150, 200],
                abundance_profiles=["lognormal", "exponential"],
                threads=1,
            )

    def test_simulate_reads_mason_expand_arguments(self):
        f1 = MagicMock(side_effect=(("reads1",), ("reads2",)))
        f2 = MagicMock(return_value=("collated_reads",))
        mock_action = MagicMock(side_effect=[f1, f2])
        mock_ctx = MagicMock(get_action=mock_action)

        simulate_reads_mason(
            ctx=mock_ctx,
            reference_genomes=self.refs,
            sample_names=["sample1", "sample2"],
            num_reads=[2000],
            read_length=[125],
            abundance_profiles=["uniform"],
        )

        f1.assert_has_calls(
            [
                call(
                    reference_genomes=self.refs,
                    sample_name="sample1",
                    num_reads=2000,
                    read_length=125,
                    abundance_profile="uniform",
                    random_seed=42,
                    threads=1,
                ),
                call(
                    reference_genomes=self.refs,
                    sample_name="sample2",
                    num_reads=2000,
                    read_length=125,
                    abundance_profile="uniform",
                    random_seed=42,
                    threads=1,
                ),
            ]
        )
        f2.assert_called_once_with(["reads1", "reads2"])


class TestGenerateAbundances(TestPluginBase):
    package = "q2_assembly.tests"

    def test_uniform_profile(self):
        result = generate_abundances("uniform", 4)
        expected = [0.25, 0.25, 0.25, 0.25]
        for r, e in zip(result, expected):
            self.assertAlmostEqual(r, e)

    def test_exponential_profile(self):
        result = generate_abundances("exponential", 3, lambd=1.0)
        # Should sum to 1
        self.assertAlmostEqual(sum(result), 1.0)
        self.assertTrue(all(x > 0 for x in result))

    def test_lognormal_profile(self):
        result = generate_abundances("lognormal", 5, mu=0, sigma=1)
        self.assertAlmostEqual(sum(result), 1.0)
        self.assertEqual(len(result), 5)

    def test_invalid_profile(self):
        with self.assertRaisesRegex(
            ValueError, "Invalid abundance profile option: 'invalid'"
        ):
            generate_abundances("invalid", 2)


class TestAbundancesToDF(TestPluginBase):
    package = "q2_assembly.tests"

    def test_abundances_to_df(self):
        abundances = [0.1, 0.2, 0.7]
        genome_files = [
            "/tmp/genomeA.fasta",
            "/tmp/genomeB.fasta",
            "/tmp/genomeC.fasta",
        ]
        sample_id = "sample1"
        df = abundances_to_df(abundances, genome_files, sample_id)
        self.assertListEqual(list(df.index), ["genomeA", "genomeB", "genomeC"])
        self.assertListEqual(list(df[sample_id]), abundances)

    def test_abundances_to_df_empty(self):
        abundances = []
        genome_files = []
        sample_id = "sample1"
        df = abundances_to_df(abundances, genome_files, sample_id)
        self.assertEqual(df.shape, (0, 1))
        self.assertListEqual(list(df.columns), [sample_id])


class TestCombineReadsAndProcessSample(TestPluginBase):
    package = "q2_assembly.tests"

    def test_combine_reads(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_id = "sampleX"
            # Create fake files to combine
            f1 = os.path.join(tmpdir, f"{sample_id}_A_00_L001_R1_001.fastq.gz")
            f2 = os.path.join(tmpdir, f"{sample_id}_B_00_L001_R1_001.fastq.gz")
            for f in [f1, f2]:
                with open(f, "wb") as fh:
                    fh.write(b"testdata\n")

            _combine_reads(sample_id, tmpdir, orientation="forward")

            # Output file should exist
            out_file = os.path.join(tmpdir, f"{sample_id}_00_L001_R1_001.fastq.gz")
            self.assertTrue(os.path.exists(out_file))

            # Output file content should be the concatenation of the two input files
            with open(out_file, "rb") as fh:
                content = fh.read()
            self.assertEqual(content, b"testdata\ntestdata\n")

            # Input files should be removed
            self.assertFalse(os.path.exists(f1))
            self.assertFalse(os.path.exists(f2))

    @patch("q2_assembly.mason.run_command")
    @patch("q2_assembly.mason._combine_reads")
    def test_process_sample(self, mock_combine_reads, mock_run_command):
        sample = "sampleY"
        genome_files = ["/tmp/genome1.fasta", "/tmp/genome2.fasta"]
        abundances = [0.6, 0.4]
        total_reads = 1000
        tmp_dir = tempfile.mkdtemp()
        threads = 2
        read_len = 150
        seed = 123

        _process_sample(
            sample,
            genome_files,
            abundances,
            total_reads,
            tmp_dir,
            threads,
            read_len,
            seed,
        )
        # Should call run_command for each genome with the expected command
        expected_calls = []
        for genome_file, abundance in zip(genome_files, abundances):
            genome_reads = int(total_reads * abundance)
            _id = os.path.basename(genome_file).replace(".fasta", "")
            expected_cmd = [
                "mason_simulator",
                "-v",
                "--seed",
                str(seed),
                "--illumina-read-length",
                str(read_len),
                "--seq-technology",
                "illumina",
                "--num-threads",
                str(threads),
                "--input-reference",
                genome_file,
                "--num-fragments",
                str(genome_reads),
                "--out",
                os.path.join(tmp_dir, f"{sample}_{_id}_00_L001_R1_001.fastq.gz"),
                "--out-right",
                os.path.join(tmp_dir, f"{sample}_{_id}_00_L001_R2_001.fastq.gz"),
            ]
            expected_calls.append(call(expected_cmd, verbose=True))

        mock_run_command.assert_has_calls(expected_calls, any_order=False)

        # Should call _combine_reads twice (forward and reverse)
        mock_combine_reads.assert_has_calls(
            [
                call(sample, tmp_dir, orientation="forward"),
                call(sample, tmp_dir, orientation="reverse"),
            ]
        )


if __name__ == "__main__":
    unittest.main()
