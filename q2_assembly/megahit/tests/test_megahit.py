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
from unittest.mock import ANY, call, patch

from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt
from qiime2 import Artifact
from qiime2.plugin.testing import TestPluginBase
from qiime2.sdk.parallel_config import ParallelConfig

from q2_assembly.megahit.megahit import (
    _assemble_megahit,
    _process_megahit_arg,
    _process_sample,
    assemble_megahit_helper,
)


class MockTempDir(tempfile.TemporaryDirectory):
    pass


class TestMegahit(TestPluginBase):
    package = "q2_assembly.tests"

    def setUp(self):
        super().setUp()
        self.assemble_megahit = self.plugin.pipelines["assemble_megahit"]
        self.fake_common_args = ["--presets", "meta-fake", "--k-min", "39"]
        self.test_params_dict = {
            "presets": "meta-sensitive",
            "min_count": 3,
            "k_list": [15, 35],
            "k_min": 17,
            "k_max": 73,
            "k_step": 14,
            "no_mercy": True,
            "bubble_level": 1,
            "prune_level": 1,
            "prune_depth": 3,
            "disconnect_ratio": 0.6,
            "low_local_ratio": 0.1,
            "max_tip_len": 3,
            "cleaning_rounds": 4,
            "no_local": False,
            "kmin_1pass": True,
            "memory": 0.75,
            "mem_flag": 1,
            "num_cpu_threads": 2,
            "no_hw_accel": False,
            "min_contig_len": 100,
        }
        self.test_params_list = [
            "--presets",
            "meta-sensitive",
            "--min-count",
            "3",
            "--k-list",
            "15,35",
            "--k-min",
            "17",
            "--k-max",
            "73",
            "--k-step",
            "14",
            "--no-mercy",
            "--bubble-level",
            "1",
            "--prune-level",
            "1",
            "--prune-depth",
            "3",
            "--disconnect-ratio",
            "0.6",
            "--low-local-ratio",
            "0.1",
            "--max-tip-len",
            "3",
            "--cleaning-rounds",
            "4",
            "--kmin-1pass",
            "--memory",
            "0.75",
            "--mem-flag",
            "1",
            "--num-cpu-threads",
            "2",
            "--min-contig-len",
            "100",
        ]

    def get_reads_path(self, kind="paired", sample_id=1, direction="fwd"):
        d = 1 if direction == "fwd" else 2
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

    def test_process_megahit_arg_simple1(self):
        obs = _process_megahit_arg("not_k_list", 123)
        exp = ["--not-k-list", "123"]
        self.assertListEqual(obs, exp)

    def test_process_megahit_arg_simple2(self):
        obs = _process_megahit_arg("k_list", [1, 2, 3])
        exp = ["--k-list", "1,2,3"]
        self.assertListEqual(obs, exp)

    def test_process_megahit_arg_bool(self):
        obs = _process_megahit_arg("k_bool", True)
        exp = ["--k-bool"]
        self.assertListEqual(obs, exp)

    @patch("subprocess.run")
    @patch("tempfile.TemporaryDirectory")
    def test_process_sample_single_end(self, p1, p2):
        result = SingleLanePerSampleSingleEndFastqDirFmt()
        contigs = self.get_data_path("sample_contigs.fa")

        test_temp_dir = MockTempDir()
        os.mkdir(os.path.join(test_temp_dir.name, "results"))
        shutil.copy(
            contigs, os.path.join(test_temp_dir.name, "results", "final.contigs.fa")
        )
        p1.return_value = test_temp_dir
        _process_sample(
            "test_sample", "fwd_reads.fastq.gz", None, self.fake_common_args, result
        )

        exp_cmd = [
            "megahit",
            "-r",
            "fwd_reads.fastq.gz",
            "-o",
            os.path.join(test_temp_dir.name, "results"),
            "--presets",
            "meta-fake",
            "--k-min",
            "39",
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
            contigs, os.path.join(test_temp_dir.name, "results", "final.contigs.fa")
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
            "megahit",
            "-1",
            "fwd_reads.fastq.gz",
            "-2",
            "rev_reads.fastq.gz",
            "-o",
            os.path.join(test_temp_dir.name, "results"),
            "--presets",
            "meta-fake",
            "--k-min",
            "39",
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
            Exception, "An error.*while running MEGAHIT.*code 123"
        ):
            _process_sample(
                "test_sample", "fwd_reads.fastq.gz", None, self.fake_common_args, result
            )

    @patch("q2_assembly.megahit._process_sample")
    def test_assemble_megahit_paired(self, p):
        input_files = self.get_data_path("reads/paired-end")
        input = SingleLanePerSamplePairedEndFastqDirFmt(input_files, mode="r")

        obs = assemble_megahit_helper(seqs=input, common_args=self.test_params_list)
        exp_calls = self.generate_exp_calls(sample_ids=(1, 2), kind="paired")

        p.assert_has_calls(exp_calls, any_order=False)
        self.assertIsInstance(obs, ContigSequencesDirFmt)

    @patch("q2_assembly.megahit._process_sample")
    def test_assemble_megahit_single(self, p):
        input_files = self.get_data_path("reads/single-end")
        input = SingleLanePerSampleSingleEndFastqDirFmt(input_files, mode="r")

        obs = assemble_megahit_helper(seqs=input, common_args=self.test_params_list)
        exp_calls = self.generate_exp_calls(sample_ids=(1, 2), kind="single")

        p.assert_has_calls(exp_calls, any_order=False)
        self.assertIsInstance(obs, ContigSequencesDirFmt)

    @patch("q2_assembly.megahit.assemble_megahit_helper")
    def test_assemble_megahit_process_params(self, p):
        input_files = self.get_data_path("reads/single-end")
        input = SingleLanePerSampleSingleEndFastqDirFmt(input_files, mode="r")

        _ = _assemble_megahit(
            seqs=input,
            presets="meta-sensitive",
            bubble_level=1,
            k_list=[1, 2],
            no_mercy=True,
        )
        exp_args = [
            "--presets",
            "meta-sensitive",
            "--min-count",
            "2",
            "--k-list",
            "1,2",
            "--no-mercy",
            "--bubble-level",
            "1",
            "--prune-level",
            "2",
            "--prune-depth",
            "2",
            "--disconnect-ratio",
            "0.1",
            "--low-local-ratio",
            "0.2",
            "--cleaning-rounds",
            "5",
            "--memory",
            "0.9",
            "--mem-flag",
            "1",
            "--num-cpu-threads",
            "1",
            "--min-contig-len",
            "200",
        ]
        p.assert_called_with(seqs=input, common_args=exp_args)

    def test_assemble_megahit_parallel_paired(self):
        input_files = self.get_data_path("formatted-reads/paired-end")
        _input = SingleLanePerSamplePairedEndFastqDirFmt(input_files, mode="r")
        samples = Artifact.import_data(
            "SampleData[PairedEndSequencesWithQuality]", _input
        )

        with ParallelConfig():
            (out,) = self.assemble_megahit.parallel(samples)._result()

        out.validate()
        self.assertIs(out.format, ContigSequencesDirFmt)

    def test_assemble_megahit_parallel_single(self):
        input_files = self.get_data_path("formatted-reads/single-end")
        _input = SingleLanePerSampleSingleEndFastqDirFmt(input_files, mode="r")
        samples = Artifact.import_data("SampleData[SequencesWithQuality]", _input)

        with ParallelConfig():
            (out,) = self.assemble_megahit.parallel(samples)._result()

        out.validate()
        self.assertIs(out.format, ContigSequencesDirFmt)


if __name__ == "__main__":
    unittest.main()
