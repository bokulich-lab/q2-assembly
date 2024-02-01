# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import shutil
import tempfile
import unittest
from subprocess import CalledProcessError
from unittest.mock import ANY, MagicMock, call, patch

from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt
from qiime2.plugin.testing import TestPluginBase

from q2_assembly.spades.spades import (
    _assemble_spades,
    _process_sample,
    _process_spades_arg,
    assemble_spades,
)


class MockTempDir(tempfile.TemporaryDirectory):
    pass


class TestSpades(TestPluginBase):
    package = "q2_assembly.tests"

    def setUp(self):
        super().setUp()
        self.fake_common_args = ["--meta", "--threads", "8"]
        self.test_params_dict = {
            "isolate": False,
            "sc": False,
            "meta": True,
            "bio": False,
            "corona": False,
            "plasmid": True,
            "metaviral": False,
            "metaplasmid": False,
            "only_assembler": True,
            "careful": True,
            "disable_rr": False,
            "threads": 8,
            "memory": 125,
            "debug": True,
            "k": [19, 27, 39],
            "cov_cutoff": 0.6,
            "phred_offset": 33,
        }
        self.test_params_list = [
            "--meta",
            "--plasmid",
            "--only-assembler",
            "--careful",
            "--threads",
            "8",
            "--memory",
            "125",
            "-k",
            "19,27,39",
            "--cov-cutoff",
            "0.6",
            "--phred-offset",
            "33",
            "--debug",
        ]

    def get_reads_path(
        self, kind="paired", sample_id=1, direction="fwd", is_single_sample=False
    ):
        d = 1 if direction == "fwd" else 2
        if is_single_sample:
            return self.get_data_path(
                f"reads/single-sample/{kind}-end/reads{sample_id}_R{d}.fastq.gz"
            )
        return self.get_data_path(f"reads/{kind}-end/reads{sample_id}_R{d}.fastq.gz")

    def generate_exp_calls(self, sample_ids, kind="paired"):
        exp_calls = []
        rev = None
        for s in sample_ids:
            fwd = self.get_reads_path(kind, s, "fwd")
            if kind == "paired":
                rev = self.get_reads_path(kind, s, "rev")
            exp_calls.append(call(f"sample{s}", fwd, rev, self.test_params_list, ANY))
        return exp_calls

    def test_process_spades_arg_simple1(self):
        obs = _process_spades_arg("not_k_list", 123)
        exp = ["--not-k-list", "123"]
        self.assertListEqual(obs, exp)

    def test_process_spades_arg_simple2(self):
        obs = _process_spades_arg("k_list", [1, 2, 3])
        exp = ["--k-list", "1,2,3"]
        self.assertListEqual(obs, exp)

    def test_process_spades_arg_k(self):
        obs = _process_spades_arg("k", [1, 2, 3])
        exp = ["-k", "1,2,3"]
        self.assertListEqual(obs, exp)

    def test_process_spades_arg_bool(self):
        obs = _process_spades_arg("k_bool", True)
        exp = ["--k-bool"]
        self.assertListEqual(obs, exp)

    def test_process_spades_arg_phred_offset_auto(self):
        obs = _process_spades_arg("phred_offset", "auto-detect")
        exp = []
        self.assertListEqual(obs, exp)

    def test_process_spades_arg_phred_offset_int(self):
        obs = _process_spades_arg("phred_offset", 33)
        exp = ["--phred-offset", "33"]
        self.assertListEqual(obs, exp)

    @patch("subprocess.run")
    @patch("tempfile.TemporaryDirectory")
    def test_process_sample_single_end(self, p1, p2):
        result = SingleLanePerSampleSingleEndFastqDirFmt()
        contigs = self.get_data_path("sample_contigs.fa")

        test_temp_dir = MockTempDir()
        os.mkdir(os.path.join(test_temp_dir.name, "results"))
        shutil.copy(
            contigs, os.path.join(test_temp_dir.name, "results", "contigs.fasta")
        )
        p1.return_value = test_temp_dir
        _process_sample(
            "test_sample", "fwd_reads.fastq.gz", None, self.fake_common_args, result
        )

        exp_cmd = [
            "spades.py",
            "-s",
            "fwd_reads.fastq.gz",
            "-o",
            os.path.join(test_temp_dir.name, "results"),
            "--meta",
            "--threads",
            "8",
        ]
        p2.assert_called_once_with(exp_cmd, check=True)

        exp_contigs = os.path.join(str(result), "test_sample_contigs.fa")
        self.assertTrue(os.path.isfile(exp_contigs))

    @patch("subprocess.run")
    @patch("tempfile.TemporaryDirectory")
    def test_process_sample_paired_end(self, p1, p2):
        result = SingleLanePerSamplePairedEndFastqDirFmt()
        contigs = self.get_data_path("sample_contigs.fa")

        test_temp_dir = MockTempDir()
        os.mkdir(os.path.join(test_temp_dir.name, "results"))
        shutil.copy(
            contigs, os.path.join(test_temp_dir.name, "results", "contigs.fasta")
        )
        p1.return_value = test_temp_dir
        _process_sample(
            "test_sample",
            "fwd_reads.fastq.gz",
            "rev_reads.fastq.gz",
            self.fake_common_args,
            result,
        )

        exp_cmd = [
            "spades.py",
            "-1",
            "fwd_reads.fastq.gz",
            "-2",
            "rev_reads.fastq.gz",
            "-o",
            os.path.join(test_temp_dir.name, "results"),
            "--meta",
            "--threads",
            "8",
        ]
        p2.assert_called_once_with(exp_cmd, check=True)

        exp_contigs = os.path.join(str(result), "test_sample_contigs.fa")
        self.assertTrue(os.path.isfile(exp_contigs))

    @patch(
        "subprocess.run", side_effect=CalledProcessError(returncode=123, cmd="some cmd")
    )
    @patch("tempfile.TemporaryDirectory")
    def test_process_sample_with_error(self, p1, p2):
        result = SingleLanePerSampleSingleEndFastqDirFmt()

        with self.assertRaisesRegex(
            Exception, "An error.*while running SPAdes.*code 123"
        ):
            _process_sample(
                "test_sample", "fwd_reads.fastq.gz", None, self.fake_common_args, result
            )

    @patch("q2_assembly.spades._process_sample")
    def test_assemble_spades_paired(self, p):
        input_files = self.get_data_path("reads/paired-end")
        input = SingleLanePerSamplePairedEndFastqDirFmt(input_files, mode="r")

        obs = _assemble_spades(
            seqs=input, meta=False, common_args=self.test_params_list
        )
        exp_calls = self.generate_exp_calls(sample_ids=(1, 2), kind="paired")

        p.assert_has_calls(exp_calls, any_order=False)
        self.assertIsInstance(obs, ContigSequencesDirFmt)

    @patch(
        "tempfile.TemporaryDirectory",
        return_value=MagicMock(__enter__=MagicMock(return_value="/tmp/mock_tmp_dir")),
    )
    @patch("q2_assembly.spades._process_sample")
    def test_assemble_spades_paired_coassemble(self, p1, p2):
        input_files = self.get_data_path("reads/paired-end")
        input = SingleLanePerSamplePairedEndFastqDirFmt(input_files, mode="r")

        obs = _assemble_spades(
            seqs=input, meta=False, coassemble=True, common_args=self.test_params_list
        )
        # to be modified accordingly
        # if reads type change in the tests/data/reads/paired_end directory
        fwd = os.path.join("/tmp", "mock_tmp_dir", "all_samples_fwd.fastq.gz")
        rev = os.path.join("/tmp", "mock_tmp_dir", "all_samples_rev.fastq.gz")
        exp_calls = [call("all_samples", fwd, rev, self.test_params_list, ANY)]

        p1.assert_has_calls(exp_calls, any_order=False)
        self.assertIsInstance(obs, ContigSequencesDirFmt)

    @patch(
        "tempfile.TemporaryDirectory",
        return_value=MagicMock(__enter__=MagicMock(return_value="/tmp/mock_tmp_dir")),
    )
    @patch("q2_assembly.spades._process_sample")
    def test_assemble_spades_paired_single_sample_coassemble(self, p1, p2):
        input_files = self.get_data_path("reads/single-sample/paired-end")
        input = SingleLanePerSamplePairedEndFastqDirFmt(input_files, mode="r")

        obs = _assemble_spades(
            seqs=input, meta=False, coassemble=True, common_args=self.test_params_list
        )
        # to be modified accordingly
        # if reads type change in the
        # tests/data/reads/single-samples/paired_end directory
        fwd = os.path.join("/tmp", "mock_tmp_dir", "all_samples_fwd.fastq.gz")
        rev = os.path.join("/tmp", "mock_tmp_dir", "all_samples_rev.fastq.gz")

        exp_calls = [call("all_samples", fwd, rev, self.test_params_list, ANY)]
        p1.assert_has_calls(exp_calls, any_order=False)
        self.assertIsInstance(obs, ContigSequencesDirFmt)

    @patch("q2_assembly.spades._process_sample")
    def test_assemble_spades_single(self, p):
        input_files = self.get_data_path("reads/single-end")
        input = SingleLanePerSampleSingleEndFastqDirFmt(input_files, mode="r")

        with self.assertRaisesRegex(
            NotImplementedError, 'SPAdes v3.15.2 in "meta" mode supports'
        ):
            _assemble_spades(seqs=input, meta=True, common_args=self.test_params_list)

    @patch("q2_assembly.spades._process_sample")
    def test_assemble_spades_single_coassemble(self, p):
        input_files = self.get_data_path("reads/single-end")
        input = SingleLanePerSampleSingleEndFastqDirFmt(input_files, mode="r")

        with self.assertRaisesRegex(
            NotImplementedError, 'SPAdes v3.15.2 in "meta" mode supports'
        ):
            _assemble_spades(
                seqs=input,
                meta=True,
                coassemble=True,
                common_args=self.test_params_list,
            )

    @patch("q2_assembly.spades._process_sample")
    def test_assemble_spades_single_single_sample_coassemble(self, p):
        input_files = self.get_data_path("reads/single-sample/single-end")
        input = SingleLanePerSampleSingleEndFastqDirFmt(input_files, mode="r")

        with self.assertRaisesRegex(
            NotImplementedError, 'SPAdes v3.15.2 in "meta" mode supports'
        ):
            _assemble_spades(
                seqs=input,
                meta=True,
                coassemble=True,
                common_args=self.test_params_list,
            )

    @patch("q2_assembly.spades._assemble_spades")
    def test_assemble_spades_process_params(self, p):
        input_files = self.get_data_path("reads/single-end")
        input = SingleLanePerSampleSingleEndFastqDirFmt(input_files, mode="r")

        _ = assemble_spades(
            seqs=input, meta=True, threads=14, k=[1, 2], cov_cutoff="off"
        )
        exp_args = [
            "--meta",
            "--threads",
            "14",
            "--memory",
            "250",
            "-k",
            "1,2",
            "--cov-cutoff",
            "off",
        ]
        p.assert_called_with(
            seqs=input, meta=True, coassemble=False, common_args=exp_args
        )


if __name__ == "__main__":
    unittest.main()
