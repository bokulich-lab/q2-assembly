# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import unittest
from unittest.mock import call, patch, MagicMock

import pandas as pd
from q2_types.genome_data import GenomeSequencesDirectoryFormat
from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt
from qiime2.plugin.testing import TestPluginBase

from q2_assembly.mason.mason import (
    _process_mason_arg,
    _simulate_reads_mason,
    generate_abundances,
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

    @patch("q2_assembly.mason.mason._process_sample")
    def test_simulate_reads_mason_helper(self, p_process):
        mock_genomes_dir_fmt = GenomeSequencesDirectoryFormat()

        refs = GenomeSequencesDirectoryFormat(
            self.get_data_path("genomes-dir-format1"), "r"
        )

        with patch(
            "q2_assembly.mason.mason.GenomeSequencesDirectoryFormat",
            return_value=mock_genomes_dir_fmt,
        ):
            reads, ft = _simulate_reads_mason(
                reference_genomes=refs,
                sample_name="sample1",
                num_reads=2000000,
                read_length=125,
                abundance_profile="uniform",
            )
        self.assertIsInstance(reads, CasavaOneEightSingleLanePerSampleDirFmt)
        p_process.assert_called_once_with(
            sample="sample1",
            genomes=mock_genomes_dir_fmt,
            abundances=pd.DataFrame(
                data={"sample1": [0.5, 0.5]},
                index=pd.Index(["ref1", "ref2"], name="id")
            ),
            total_reads=2000000,
            results_dir=str(reads),
            threads=1,
            read_len=125,
            seed=42,
        )
        expected_ft = pd.DataFrame(
            data={"ref1": [0.5], "ref2": [0.5]},
            index=pd.Index(["sample1"], name="id")
        )
        pd.testing.assert_frame_equal(ft, expected_ft)

    @patch("q2_assembly.mason.mason._process_sample")
    def test_simulate_reads_mason_helper_with_abundances(self, p_process):
        mock_genomes_dir_fmt = GenomeSequencesDirectoryFormat()

        refs = GenomeSequencesDirectoryFormat(
            self.get_data_path("genomes-dir-format1"), "r"
        )

        # Pre-calculated abundances
        abundances = pd.DataFrame(
            data={"sample1": [0.7, 0.3]},
            index=pd.Index(["ref1", "ref2"], name="id")
        )

        with patch(
            "q2_assembly.mason.mason.GenomeSequencesDirectoryFormat",
            return_value=mock_genomes_dir_fmt,
        ):
            reads, ft = _simulate_reads_mason(
                reference_genomes=refs,
                sample_name="sample1",
                num_reads=2000000,
                read_length=125,
                abundances=abundances,
            )
        self.assertIsInstance(reads, CasavaOneEightSingleLanePerSampleDirFmt)
        p_process.assert_called_once_with(
            sample="sample1",
            genomes=mock_genomes_dir_fmt,
            abundances=abundances,
            total_reads=2000000,
            results_dir=str(reads),
            threads=1,
            read_len=125,
            seed=42,
        )
        expected_ft = pd.DataFrame(
            data={"ref1": [0.7], "ref2": [0.3]},
            index=pd.Index(["sample1"], name="id")
        )
        pd.testing.assert_frame_equal(ft, expected_ft)

    def test_simulate_reads_mason_helper_mutual_exclusivity_both(self):
        refs = GenomeSequencesDirectoryFormat(
            self.get_data_path("genomes-dir-format1"), "r"
        )
        abundances = pd.DataFrame(
            data={"sample1": [0.7, 0.3]},
            index=pd.Index(["ref1", "ref2"], name="id")
        )

        with self.assertRaisesRegex(
            ValueError,
            "Cannot provide both 'abundance_profile' and 'abundances'"
        ):
            _simulate_reads_mason(
                reference_genomes=refs,
                sample_name="sample1",
                abundance_profile="uniform",
                abundances=abundances,
            )

    def test_simulate_reads_mason_helper_mutual_exclusivity_neither(self):
        refs = GenomeSequencesDirectoryFormat(
            self.get_data_path("genomes-dir-format1"), "r"
        )

        with self.assertRaisesRegex(
            ValueError,
            "Must provide either 'abundance_profile' or 'abundances'"
        ):
            _simulate_reads_mason(
                reference_genomes=refs,
                sample_name="sample1",
            )

    def test_simulate_reads_mason_one_sample(self):
        f1 = MagicMock(return_value=("reads", "ft"))
        f2 = MagicMock(return_value=("collated_reads",))
        f3 = MagicMock(return_value=("merged_tables",))
        mock_action = MagicMock(side_effect=[f1, f2, f3])
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
        f1 = MagicMock(side_effect=(("reads1", "ft1"), ("reads2", "ft2")))
        f2 = MagicMock(return_value=("collated_reads",))
        f3 = MagicMock(return_value=("merged_tables",))
        mock_action = MagicMock(side_effect=[f1, f2, f3])
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

    def test_simulate_reads_mason_with_abundances_single_sample(self):
        f1 = MagicMock(return_value=("reads", "ft"))
        f2 = MagicMock(return_value=("collated_reads",))
        f3 = MagicMock(return_value=("merged_tables",))
        mock_action = MagicMock(side_effect=[f1, f2, f3])
        mock_ctx = MagicMock(get_action=mock_action)

        abundances = pd.DataFrame(
            data={"sample1": [0.6, 0.4]},
            index=pd.Index(["ref1", "ref2"], name="id")
        )

        simulate_reads_mason(
            ctx=mock_ctx,
            reference_genomes=self.refs,
            abundances=abundances,
            num_reads=[2000000],
            read_length=[125],
        )

        # Verify that abundances for sample1 were passed
        f1.assert_called_once_with(
            reference_genomes=self.refs,
            sample_name="sample1",
            num_reads=2000000,
            read_length=125,
            abundances=abundances[["sample1"]],
            random_seed=42,
            threads=1,
        )
        f2.assert_called_once_with(["reads"])

    def test_simulate_reads_mason_with_abundances_multiple_samples(self):
        f1 = MagicMock(side_effect=(("reads1", "ft1"), ("reads2", "ft2")))
        f2 = MagicMock(return_value=("collated_reads",))
        f3 = MagicMock(return_value=("merged_tables",))
        mock_action = MagicMock(side_effect=[f1, f2, f3])
        mock_ctx = MagicMock(get_action=mock_action)

        abundances = pd.DataFrame(
            data={"sample1": [0.6, 0.4], "sample2": [0.3, 0.7]},
            index=pd.Index(["ref1", "ref2"], name="id")
        )

        simulate_reads_mason(
            ctx=mock_ctx,
            reference_genomes=self.refs,
            abundances=abundances,
            num_reads=[2000, 4000],
            read_length=[125, 150],
        )

        f1.assert_has_calls(
            [
                call(
                    reference_genomes=self.refs,
                    sample_name="sample1",
                    num_reads=2000,
                    read_length=125,
                    abundances=abundances[["sample1"]],
                    random_seed=42,
                    threads=1,
                ),
                call(
                    reference_genomes=self.refs,
                    sample_name="sample2",
                    num_reads=4000,
                    read_length=150,
                    abundances=abundances[["sample2"]],
                    random_seed=42,
                    threads=1,
                ),
            ]
        )
        f2.assert_called_once_with(["reads1", "reads2"])

    def test_simulate_reads_mason_mutual_exclusivity_both(self):
        abundances = pd.DataFrame(
            data={"sample1": [0.6, 0.4]},
            index=pd.Index(["ref1", "ref2"], name="id")
        )

        with self.assertRaisesRegex(
            ValueError,
            "Cannot provide both 'abundance_profiles' and 'abundances'"
        ):
            simulate_reads_mason(
                ctx=MagicMock(),
                reference_genomes=self.refs,
                abundance_profiles=["uniform"],
                abundances=abundances,
            )

    def test_simulate_reads_mason_mutual_exclusivity_neither(self):
        with self.assertRaisesRegex(
            ValueError,
            "Must provide either 'abundance_profiles' or 'abundances'"
        ):
            simulate_reads_mason(
                ctx=MagicMock(),
                reference_genomes=self.refs,
                sample_names=["sample1"],
            )

    def test_simulate_reads_mason_sample_names_with_abundances_error(self):
        abundances = pd.DataFrame(
            data={"sample1": [0.6, 0.4]},
            index=pd.Index(["ref1", "ref2"], name="id")
        )

        with self.assertRaisesRegex(
            ValueError,
            "Cannot provide 'sample_names' when 'abundances' is provided"
        ):
            simulate_reads_mason(
                ctx=MagicMock(),
                reference_genomes=self.refs,
                sample_names=["sample1"],
                abundances=abundances,
            )

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
        f1 = MagicMock(side_effect=(("reads1", "ft1"), ("reads2", "ft2")))
        f2 = MagicMock(return_value=("collated_reads",))
        f3 = MagicMock(return_value=("merged_tables",))
        mock_action = MagicMock(side_effect=[f1, f2, f3])
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

    def test_simulate_reads_mason_abundances_extract_sample_names(self):
        f1 = MagicMock(side_effect=(("reads1", "ft1"), ("reads2", "ft2"), ("reads3", "ft3")))
        f2 = MagicMock(return_value=("collated_reads",))
        f3 = MagicMock(return_value=("merged_tables",))
        mock_action = MagicMock(side_effect=[f1, f2, f3])
        mock_ctx = MagicMock(get_action=mock_action)

        # Abundances with three samples
        abundances = pd.DataFrame(
            data={
                "sampleA": [0.5, 0.5],
                "sampleB": [0.3, 0.7],
                "sampleC": [0.8, 0.2],
            },
            index=pd.Index(["ref1", "ref2"], name="id")
        )

        simulate_reads_mason(
            ctx=mock_ctx,
            reference_genomes=self.refs,
            abundances=abundances,
            num_reads=[1000],
            read_length=[100],
        )

        # Should extract sample names from columns
        f1.assert_has_calls(
            [
                call(
                    reference_genomes=self.refs,
                    sample_name="sampleA",
                    num_reads=1000,
                    read_length=100,
                    abundances=abundances[["sampleA"]],
                    random_seed=42,
                    threads=1,
                ),
                call(
                    reference_genomes=self.refs,
                    sample_name="sampleB",
                    num_reads=1000,
                    read_length=100,
                    abundances=abundances[["sampleB"]],
                    random_seed=42,
                    threads=1,
                ),
                call(
                    reference_genomes=self.refs,
                    sample_name="sampleC",
                    num_reads=1000,
                    read_length=100,
                    abundances=abundances[["sampleC"]],
                    random_seed=42,
                    threads=1,
                ),
            ]
        )
        f2.assert_called_once_with(["reads1", "reads2", "reads3"])


class TestGenerateAbundances(TestPluginBase):
    package = "q2_assembly.tests"

    def test_uniform_profile(self):
        genomes = GenomeSequencesDirectoryFormat(self.get_data_path("genomes-dir-format4"), "r")
        result = generate_abundances("uniform", genomes, "sample1")
        expected = pd.DataFrame(
            data={"sample1": [0.25, 0.25, 0.25, 0.25]},
            index=pd.Index(["ref1", "ref2", "ref3", "ref4"], name="id")
        )
        pd.testing.assert_frame_equal(result, expected)

    def test_exponential_profile(self):
        genomes = GenomeSequencesDirectoryFormat(self.get_data_path("genomes-dir-format3"), "r")
        result = generate_abundances("exponential", genomes, "sample1", lambd=1.0)
        expected = pd.DataFrame(
            data={"sample1": [0.6652, 0.2447, 0.0900]},
            index=pd.Index(["ref1", "ref2", "ref3"], name="id")
        )
        pd.testing.assert_frame_equal(result, expected, check_exact=False, rtol=0.001)

    def test_lognormal_profile(self):
        genomes = GenomeSequencesDirectoryFormat(self.get_data_path("genomes-dir-format5"), "r")
        result = generate_abundances("lognormal", genomes, "sample1", mu=0, sigma=1)
        expected = pd.DataFrame(
            data={"sample1": [0.1676, 0.0888, 0.1950, 0.4678, 0.0807]},
            index=pd.Index(["ref1", "ref2", "ref3", "ref4", "ref5"], name="id")
        )
        pd.testing.assert_frame_equal(result, expected, check_exact=False, rtol=0.001)


    def test_invalid_profile(self):
        genomes = GenomeSequencesDirectoryFormat(self.get_data_path("genomes-dir-format4"), "r")
        with self.assertRaisesRegex(
            ValueError, "Invalid abundance profile option: 'invalid'"
        ):
            generate_abundances("invalid", genomes, "sample1")


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
        genomes = GenomeSequencesDirectoryFormat(self.get_data_path("genomes-dir-format1"), "r")
        abundances = pd.DataFrame(
            data={sample: [0.6, 0.4]},
            index=pd.Index(["ref1", "ref2"], name="id")
        )
        total_reads = 1000
        tmp_dir = tempfile.mkdtemp()
        threads = 2
        read_len = 150
        seed = 123

        _process_sample(
            sample,
            genomes,
            abundances,
            total_reads,
            tmp_dir,
            threads,
            read_len,
            seed,
        )
        # Should call run_command for each genome with the expected command
        expected_calls = []
        for genome_id, genome_fp in genomes.file_dict().items():
            abundance = abundances.loc[genome_id, sample]
            genome_reads = int(total_reads * abundance)
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
                genome_fp,
                "--num-fragments",
                str(genome_reads),
                "--out",
                os.path.join(tmp_dir, f"{sample}_{genome_id}_00_L001_R1_001.fastq.gz"),
                "--out-right",
                os.path.join(tmp_dir, f"{sample}_{genome_id}_00_L001_R2_001.fastq.gz"),
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
