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
from unittest.mock import call, patch, ANY

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
)


class MockTempDir(tempfile.TemporaryDirectory):
    pass


class TestMason(TestPluginBase):
    package = "q2_assembly.tests"

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

    @patch("tempfile.TemporaryDirectory")
    @patch("q2_assembly.mason._process_sample")
    def test_simulate_reads_mason(self, p_process, p_tmpdir):
        tmp_dir = tempfile.mkdtemp()
        p_tmpdir.return_value.__enter__.return_value = tmp_dir
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
                sample_names=["sample1", "sample2"],
                num_reads=[1000000, 2000000],
                read_length=[100, 150],
                abundance_profiles=["uniform", "lognormal"],
            )
        self.assertIsInstance(result, CasavaOneEightSingleLanePerSampleDirFmt)
        p_process.assert_has_calls(
            [
                call(
                    "sample1",
                    [
                        os.path.join(str(mock_genomes_dir_fmt), "ref1.fasta"),
                        os.path.join(str(mock_genomes_dir_fmt), "ref2.fasta"),
                    ],
                    [0.5, 0.5],
                    1000000,
                    tmp_dir,
                    1,
                    100,
                    42,
                ),
                call(
                    "sample2",
                    [
                        os.path.join(str(mock_genomes_dir_fmt), "ref1.fasta"),
                        os.path.join(str(mock_genomes_dir_fmt), "ref2.fasta"),
                    ],
                    [0.6536174529063914, 0.3463825470936087],
                    2000000,
                    tmp_dir,
                    1,
                    150,
                    42,
                ),
            ]
        )

    def test_simulate_reads_mason_duplicate_sample_names(self):
        refs = GenomeSequencesDirectoryFormat(
            self.get_data_path("genomes-dir-format1"), "r"
        )
        with self.assertRaisesRegex(ValueError, "Sample names need to be unique"):
            _simulate_reads_mason(
                reference_genomes=refs,
                sample_names=["sample1", "sample1"],
                num_reads=[1000000, 2000000],
                read_length=[100, 150],
                abundance_profiles=["uniform", "lognormal"],
                threads=1,
            )

    def test_simulate_reads_mason_num_reads_length(self):
        refs = GenomeSequencesDirectoryFormat(
            self.get_data_path("genomes-dir-format1"), "r"
        )
        with self.assertRaisesRegex(
            ValueError,
            "The length of read_length must be either 1 "
            "or equal to the number of sample names",
        ):
            _simulate_reads_mason(
                reference_genomes=refs,
                sample_names=["sample1", "sample2"],
                num_reads=[1000000],
                read_length=[100, 150, 200],
                abundance_profiles=["uniform", "lognormal"],
                threads=1,
            )

    def test_simulate_reads_mason_read_length_length(self):
        refs = GenomeSequencesDirectoryFormat(
            self.get_data_path("genomes-dir-format1"), "r"
        )
        with self.assertRaisesRegex(
            ValueError,
            "The length of num_reads must be either 1 "
            "or equal to the number of sample names",
        ):
            _simulate_reads_mason(
                reference_genomes=refs,
                sample_names=["sample1", "sample2", "sample3"],
                num_reads=[1000000, 2000000],
                read_length=[100, 150, 200],
                abundance_profiles=["uniform", "lognormal", "exponential"],
                threads=1,
            )

    @patch("q2_assembly.mason._process_sample")
    @patch("tempfile.TemporaryDirectory")
    def test_simulate_reads_mason_expand_num_reads_and_read_length(
        self, p_tmpdir, p_process
    ):
        tmp_dir = tempfile.mkdtemp()
        p_tmpdir.return_value.__enter__.return_value = tmp_dir
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
                sample_names=["sample1", "sample2"],
                num_reads=[1000000],
                read_length=[100],
                abundance_profiles=["uniform", "lognormal"],
                threads=1,
            )
        self.assertIsInstance(result, CasavaOneEightSingleLanePerSampleDirFmt)
        calls = p_process.call_args_list
        self.assertEqual(len(calls), 2)
        self.assertEqual(calls[0][0][3], 1000000)
        self.assertEqual(calls[1][0][3], 1000000)
        self.assertEqual(calls[0][0][6], 100)
        self.assertEqual(calls[1][0][6], 100)


class TestGenerateAbundances(TestPluginBase):
    package = "q2_assembly.tests"

    def test_uniform_profile(self):
        profiles = ["uniform"]
        num_genomes = 4
        result = generate_abundances(profiles, num_genomes)
        expected = [[0.25, 0.25, 0.25, 0.25]]
        for r, e in zip(result[0], expected[0]):
            self.assertAlmostEqual(r, e)

    def test_exponential_profile(self):
        profiles = ["exponential"]
        num_genomes = 3
        result = generate_abundances(profiles, num_genomes, lambd=1.0)
        # Should sum to 1
        self.assertAlmostEqual(sum(result[0]), 1.0)
        self.assertTrue(all(x > 0 for x in result[0]))

    def test_lognormal_profile(self):
        profiles = ["lognormal"]
        num_genomes = 5
        result = generate_abundances(profiles, num_genomes, mu=0, sigma=1)
        self.assertAlmostEqual(sum(result[0]), 1.0)
        self.assertEqual(len(result[0]), num_genomes)

    def test_invalid_profile(self):
        profiles = ["invalid"]
        num_genomes = 2
        result = generate_abundances(profiles, num_genomes)
        self.assertListEqual(result, [[]])

    def test_multiple_profiles(self):
        profiles = ["uniform", "exponential", "lognormal"]
        num_genomes = 4
        result = generate_abundances(profiles, num_genomes)
        self.assertEqual(len(result), 3)
        for abundances in result:
            self.assertEqual(len(abundances), num_genomes)
            self.assertAlmostEqual(sum(abundances), 1.0)


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

    @patch("subprocess.run")
    def test_combine_reads(self, mock_run):
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

            # Input files should be removed
            self.assertFalse(os.path.exists(f1))
            self.assertFalse(os.path.exists(f2))

            # subprocess.run should have been called with the expected command
            expected_cmd = ["cat", f1, f2]
            mock_run.assert_called_with(expected_cmd, check=True, stdout=ANY)

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
