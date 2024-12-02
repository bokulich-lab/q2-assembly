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
from subprocess import CalledProcessError
from unittest.mock import call, patch

from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt

from q2_assembly.mason.mason import (
    _process_mason_arg,
    _simulate_reads_mason,
    simulate_reads_mason,
)


class MockTempDir(tempfile.TemporaryDirectory):
    pass


class TestMason(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.fake_common_args = ["--num_reads", "1000000", "--read_length", "100"]
        self.test_params_dict = {
            "reference_genomes": "test_genome.fa",
            "sample_names": ["sample1", "sample2"],
            "num_reads": [1000000, 2000000],
            "read_length": [100, 150],
            "fragment_mean_size": 500,
            "fragment_size_stddev": 50,
            "error_rate": 0.01,
            "random_seed": 42,
            "abundance_profiles": ["uniform", "lognormal"],
            "threads": 1,
        }
        self.test_params_list = [
            "--fragment_mean_size",
            "500",
            "--fragment_size_stddev",
            "50",
            "--error_rate",
            "0.01",
            "--random_seed",
            "42",
        ]

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

    @patch("subprocess.run")
    def test_simulate_reads_mason(self, p):
        samples = ["sample1", "sample2"]
        result_fp = "/there"
        _simulate_reads_mason(
            reference_genomes="test_genome.fa",
            sample_names=samples,
            num_reads=[1000000, 2000000],
            read_length=[100, 150],
            fragment_mean_size=500,
            fragment_size_stddev=50,
            error_rate=0.01,
            random_seed=42,
            abundance_profiles=["uniform", "lognormal"],
            threads=1,
        )

        exp_calls = []
        for s in samples:
            call_ext = self.test_params_list.copy()
            exp_calls.append(call(call_ext, check=True))
        p.assert_has_calls(exp_calls)

    @patch(
        "subprocess.run", side_effect=CalledProcessError(returncode=123, cmd="some cmd")
    )
    def test_simulate_reads_mason_with_error(self, p):
        samples = ["sample1", "sample2"]
        result_fp = "/there"

        with self.assertRaisesRegex(
            Exception, r"An error.*while running Mason.*code 123"
        ):
            _simulate_reads_mason(
                reference_genomes="test_genome.fa",
                sample_names=samples,
                num_reads=[1000000, 2000000],
                read_length=[100, 150],
                fragment_mean_size=500,
                fragment_size_stddev=50,
                error_rate=0.01,
                random_seed=42,
                abundance_profiles=["uniform", "lognormal"],
                threads=1,
            )

    @patch("shutil.move")
    @patch("q2_assembly.mason.mason._simulate_reads_mason")
    @patch("tempfile.TemporaryDirectory")
    def test_simulate_reads_mason_action(self, p1, p2, p3):
        test_temp_dir = MockTempDir()
        p1.return_value = test_temp_dir
        for x in ["sample1_00_L001_R1_001.fastq.gz", "sample1_00_L001_R2_001.fastq.gz"]:
            with open(os.path.join(test_temp_dir.name, x), "w") as f:
                f.write("test")

        obs_reads = simulate_reads_mason(**self.test_params_dict)
        p2.assert_called_with(
            reference_genomes="test_genome.fa",
            sample_names=["sample1", "sample2"],
            num_reads=[1000000, 2000000],
            read_length=[100, 150],
            fragment_mean_size=500,
            fragment_size_stddev=50,
            error_rate=0.01,
            random_seed=42,
            abundance_profiles=["uniform", "lognormal"],
            threads=1,
        )

        self.assertIsInstance(obs_reads, CasavaOneEightSingleLanePerSampleDirFmt)
        self.assertTrue(os.path.exists(str(obs_reads)))


if __name__ == "__main__":
    unittest.main()
