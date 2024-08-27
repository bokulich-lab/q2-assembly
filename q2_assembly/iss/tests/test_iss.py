# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import os
import shutil
import tempfile
import unittest
from copy import deepcopy
from subprocess import CalledProcessError
from unittest.mock import call, patch

import biom
import pandas as pd
from q2_types.feature_data import DNAIterator
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin.util import transform

from q2_assembly.iss.iss import (
    _abundances_to_biom,
    _ensure_sample_names_exists,
    _generate_reads,
    _process_iss_arg,
    _rename_reads_files,
    generate_reads,
)


class MockTempDir(tempfile.TemporaryDirectory):
    pass


class TestISS(TestPluginBase):
    package = "q2_assembly.iss.tests"

    def setUp(self):
        super().setUp()
        self.fake_common_args = ["--presets", "meta-fake", "--k-min", "39"]
        self.test_params_dict = {
            "n_genomes": 10,
            "ncbi": ["bacteria", "archaea"],
            "n_genomes_ncbi": [3, 2],
            "abundance": "halfnormal",
            "debug": True,
        }
        self.test_params_list = [
            "--n_genomes",
            "10",
            "--ncbi",
            "bacteria",
            "archaea",
            "--n_genomes_ncbi",
            "3",
            "2",
            "--abundance",
            "halfnormal",
            "--n_reads",
            "1000000",
            "--mode",
            "kde",
            "--model",
            "HiSeq",
            "--cpus",
            "1",
            "--debug",
        ]

    def generate_exp_calls(self, samples, dest):
        base_call = ["iss", "generate", "--compress"]
        base_call.extend(self.test_params_list)
        calls = []
        for s in samples:
            call_ext = deepcopy(base_call)
            call_ext.extend(["--output", os.path.join(dest, f"{s}_00_L001")])
            calls.append(call(call_ext, check=True))
        return calls

    def copy_test_files(self, dp, extension, dest):
        for f in glob.glob(os.path.join(self.get_data_path(dp), extension)):
            shutil.copy(f, dest)

    def test_process_iss_arg(self):
        obs = _process_iss_arg("not_k_list", 123)
        exp = ["--not_k_list", "123"]
        self.assertListEqual(obs, exp)

    def test_process_iss_arg_list(self):
        obs = _process_iss_arg("k_list", [1, 2, 3])
        exp = ["--k_list", "1", "2", "3"]
        self.assertListEqual(obs, exp)

    def test_process_iss_arg_bool(self):
        obs = _process_iss_arg("k_bool", True)
        exp = ["--k_bool"]
        self.assertListEqual(obs, exp)

    @patch("builtins.print")
    def test_sample_names_empty(self, mock_print):
        result = _ensure_sample_names_exists([])

        self.assertListEqual(result, ["sample"])

        mock_print.assert_called_once_with(
            'The "--p-sample-names" option was not provided. '
            'Only one sample will be created with the prefix "sample".'
            "\n"
        )

    def test_sample_names_not_empty(self):
        result = _ensure_sample_names_exists(["sample1", "sample2"])
        self.assertListEqual(result, ["sample1", "sample2"])

    @patch("os.rename")
    def test_rename_reads(self, p):
        fp = self.get_data_path("reads")
        obs_reads = _rename_reads_files(fp)
        exp_reads = [
            os.path.join(fp, f"sample{x+1}_00_L001_R{y+1}_001_001.fastq.gz")
            for x in range(2)
            for y in range(2)
        ]
        self.assertListEqual(obs_reads, exp_reads)

    @patch("subprocess.run")
    @patch("q2_assembly.iss._rename_reads_files")
    def test_generate_reads(self, p1, p2):
        samples = ["samp1", "samp2"]
        result_fp = "/there"
        _generate_reads(samples, self.test_params_list, result_fp)

        p1.assert_called_with(result_fp)
        p2.assert_has_calls(self.generate_exp_calls(samples, result_fp))

    @patch(
        "subprocess.run", side_effect=CalledProcessError(returncode=123, cmd="some cmd")
    )
    def test_generate_reads_with_error(self, p):
        samples = ["samp1", "samp2"]
        result_fp = "/there"

        with self.assertRaisesRegex(
            Exception, r"An error.*while running InSilicoSeq.*code 123"
        ):
            _generate_reads(samples, self.test_params_list, result_fp)

    def test_abundances_to_biom(self):
        abunds_fp = sorted(
            glob.glob(os.path.join(self.get_data_path("abundances"), "*.txt"))
        )
        obs_table = _abundances_to_biom(abunds_fp)
        obs_df = transform(obs_table, to_type=pd.DataFrame)

        with open(self.get_data_path("abundances/biom_table.tsv"), "r") as f:
            exp_table = biom.Table.from_tsv(f, None, None, lambda x: x)
            exp_df = transform(exp_table, to_type=pd.DataFrame)

        pd.testing.assert_frame_equal(obs_df, exp_df)

    def test_generate_reads_wrong_genome_counts(self):
        with self.assertRaisesRegex(Exception, r".*provided 2 kingdom\(s\) but 1.*"):
            generate_reads(ncbi=["bacteria", "archaea"], n_genomes_ncbi=[5])

    def test_generate_reads_duplicated_samples(self):
        with self.assertRaisesRegex(Exception, r".*duplicated names\: s1\, s2.*"):
            generate_reads(sample_names=["s1", "s2", "s3", "s1", "s2"])

    @patch("shutil.move")
    @patch("q2_assembly.iss._generate_reads")
    @patch("tempfile.TemporaryDirectory")
    def test_generate_reads_action(self, p1, p2, p3):
        test_temp_dir = MockTempDir()
        p1.return_value = test_temp_dir
        for x, y in [
            ("reads", "*.gz"),
            ("genomes", "*.fasta"),
            ("abundances", "*.txt"),
        ]:
            self.copy_test_files(x, y, str(test_temp_dir.name))

        _, obs_genomes, obs_abundances = generate_reads(
            sample_names=["samp1", "samp2"], **self.test_params_dict
        )
        p2.assert_called_with(
            ["samp1", "samp2"], self.test_params_list, test_temp_dir.name
        )

        obs_abundances_df = transform(obs_abundances, to_type=pd.DataFrame)

        exp_genomes = {
            "genome1": "ATGCATGC",
            "genome2": "GATCGCATGA",
            "genome3": "ATGATCGCATGAGC",
            "genome4": "GCATGAGCCATGA",
        }
        for seq in obs_genomes.view(DNAIterator):
            _id = seq.metadata["id"]
            self.assertEqual(str(seq), exp_genomes[_id])

        with open(self.get_data_path("abundances/biom_table.tsv"), "r") as f:
            exp_biom_table = biom.Table.from_tsv(f, None, None, lambda x: x)
            exp_df = transform(exp_biom_table, to_type=pd.DataFrame)
        pd.testing.assert_frame_equal(obs_abundances_df, exp_df)


if __name__ == "__main__":
    unittest.main()
