# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
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

from parameterized import parameterized
from q2_types.per_sample_sequences import (
    ContigSequencesDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from qiime2 import Artifact
from qiime2.plugin.testing import TestPluginBase

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

    def get_reads_path(
        self, kind="paired", sample_id=1, direction="fwd", is_single_sample=False
    ):
        d = 1 if direction == "fwd" else 2
        if is_single_sample:
            return self.get_data_path(
                f"reads/single-sample/{kind}-end/reads{sample_id}_R{d}.fastq.gz"
            )
        return self.get_data_path(f"reads/{kind}-end/reads{sample_id}_R{d}.fastq.gz")

    def get_fwd_rev_paths(self, kind, sample_id, is_single_sample=False):
        rev = None
        fwd = self.get_reads_path(kind, sample_id, "fwd", is_single_sample)
        if kind == "paired":
            rev = self.get_reads_path(kind, sample_id, "rev", is_single_sample)

        return fwd, rev

    def generate_exp_calls_coassembly(
        self,
        sample_ids,
        kind="paired",
        coassemble=False,
        uuid_type="shortuuid",
        is_single_sample=False,
    ):
        exp_calls = []
        fwd = []
        rev = []
        if coassemble:
            for s in sample_ids:
                # collect paths in 2 lists, then make a single call
                fwd.append(self.get_fwd_rev_paths(kind, s, is_single_sample)[0])
                if kind == "paired":
                    rev.append(self.get_fwd_rev_paths(kind, s, is_single_sample)[1])
            exp_calls.append(
                call(
                    "all_contigs",
                    ",".join(fwd),
                    (",".join(rev) if len(rev) != 0 else None),
                    self.test_params_list,
                    ANY,
                )
            )
        else:
            for s in sample_ids:
                # make a list of calls, each call has separate
                # samples and respective paths
                fwd, rev = self.get_fwd_rev_paths(kind, s, is_single_sample)
                exp_calls.append(
                    call(f"sample{s}", fwd, rev, self.test_params_list, ANY)
                )

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

    @patch("q2_assembly.megahit.megahit.modify_contig_ids")
    @patch("q2_assembly.megahit.megahit._process_sample")
    def test_assemble_megahit_paired(self, p1, p2):
        input_files = self.get_data_path("reads/paired-end")
        input = SingleLanePerSamplePairedEndFastqDirFmt(input_files, mode="r")

        obs = assemble_megahit_helper(
            reads=input,
            coassemble=False,
            uuid_type="shortuuid",
            common_args=self.test_params_list,
        )
        exp_calls = self.generate_exp_calls_coassembly(
            sample_ids=(1, 2), kind="paired", coassemble=False
        )

        p1.assert_has_calls(exp_calls, any_order=False)
        p2.assert_has_calls(
            [call(os.path.join(str(obs), f"sample1_contigs.fa"), "sample1", "shortuuid"),
             call(os.path.join(str(obs), f"sample2_contigs.fa"), "sample2", "shortuuid")]
        )
        self.assertIsInstance(obs, ContigSequencesDirFmt)

    @patch("q2_assembly.megahit.megahit.modify_contig_ids")
    @patch("q2_assembly.megahit.megahit._process_sample")
    def test_assemble_megahit_single(self, p1, p2):
        input_files = self.get_data_path("reads/single-end")
        input = SingleLanePerSampleSingleEndFastqDirFmt(input_files, mode="r")

        obs = assemble_megahit_helper(
            reads=input,
            coassemble=False,
            uuid_type="shortuuid",
            common_args=self.test_params_list,
        )
        exp_calls = self.generate_exp_calls_coassembly(
            sample_ids=(1, 2), kind="single", coassemble=False
        )

        p1.assert_has_calls(exp_calls, any_order=False)
        p2.assert_has_calls(
            [call(os.path.join(str(obs), f"sample1_contigs.fa"), "sample1", "shortuuid"),
             call(os.path.join(str(obs), f"sample2_contigs.fa"), "sample2", "shortuuid")]
        )
        self.assertIsInstance(obs, ContigSequencesDirFmt)

    @patch("q2_assembly.megahit.megahit.modify_contig_ids")
    @patch("q2_assembly.megahit.megahit._process_sample")
    def test_assemble_megahit_paired_coassemble(self, p1, p2):
        input_files = self.get_data_path("reads/paired-end")
        input = SingleLanePerSamplePairedEndFastqDirFmt(input_files, mode="r")

        obs = assemble_megahit_helper(
            reads=input,
            coassemble=True,
            uuid_type="shortuuid",
            common_args=self.test_params_list,
        )
        exp_calls = self.generate_exp_calls_coassembly(
            sample_ids=(1, 2), kind="paired", coassemble=True, uuid_type="shortuuid"
        )

        p1.assert_has_calls(exp_calls, any_order=False)
        p2.assert_has_calls([call(os.path.join(str(obs), f"all_contigs.fa"), "all_contigs", "shortuuid")])
        self.assertIsInstance(obs, ContigSequencesDirFmt)

    @patch("q2_assembly.megahit.megahit.modify_contig_ids")
    @patch("q2_assembly.megahit.megahit._process_sample")
    def test_assemble_megahit_single_coassemble(self, p1, p2):
        input_files = self.get_data_path("reads/single-end")
        input = SingleLanePerSampleSingleEndFastqDirFmt(input_files, mode="r")

        obs = assemble_megahit_helper(
            reads=input,
            coassemble=True,
            uuid_type="shortuuid",
            common_args=self.test_params_list,
        )
        exp_calls = self.generate_exp_calls_coassembly(
            sample_ids=(1, 2), kind="single", coassemble=True
        )

        p1.assert_has_calls(exp_calls, any_order=False)
        p2.assert_has_calls([call(os.path.join(str(obs), f"all_contigs.fa"), "all_contigs", "shortuuid")])
        self.assertIsInstance(obs, ContigSequencesDirFmt)

    @patch("q2_assembly.megahit.megahit.modify_contig_ids")
    @patch("q2_assembly.megahit.megahit._process_sample")
    def test_assemble_megahit_paired_single_sample_coassemble(self, p1, p2):
        input_files = self.get_data_path("reads/single-sample/paired-end")
        input = SingleLanePerSamplePairedEndFastqDirFmt(input_files, mode="r")

        obs = assemble_megahit_helper(
            reads=input,
            coassemble=True,
            uuid_type="shortuuid",
            common_args=self.test_params_list,
        )
        exp_calls = self.generate_exp_calls_coassembly(
            sample_ids=(1,), kind="paired", coassemble=True, is_single_sample=True
        )

        p1.assert_has_calls(exp_calls, any_order=False)
        p2.assert_has_calls([call(os.path.join(str(obs), f"all_contigs.fa"), "all_contigs", "shortuuid")])
        self.assertIsInstance(obs, ContigSequencesDirFmt)

    @patch("q2_assembly.megahit.megahit.modify_contig_ids")
    @patch("q2_assembly.megahit.megahit._process_sample")
    def test_assemble_megahit_single_single_sample_coassemble(self, p1, p2):
        input_files = self.get_data_path("reads/single-sample/single-end")
        input = SingleLanePerSampleSingleEndFastqDirFmt(input_files, mode="r")

        obs = assemble_megahit_helper(
            reads=input,
            coassemble=True,
            uuid_type="shortuuid",
            common_args=self.test_params_list,
        )
        exp_calls = self.generate_exp_calls_coassembly(
            sample_ids=(1,), kind="single", coassemble=True, is_single_sample=True
        )

        p1.assert_has_calls(exp_calls, any_order=False)
        p2.assert_has_calls([call(os.path.join(str(obs), f"all_contigs.fa"), "all_contigs", "shortuuid")])
        self.assertIsInstance(obs, ContigSequencesDirFmt)

    @patch("q2_assembly.megahit.megahit.assemble_megahit_helper")
    def test_assemble_megahit_process_params(self, p1):
        input_files = self.get_data_path("reads/single-end")
        input = SingleLanePerSampleSingleEndFastqDirFmt(input_files, mode="r")

        _ = _assemble_megahit(
            reads=input,
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
        p1.assert_called_with(
            reads=input, coassemble=False, uuid_type="shortuuid", common_args=exp_args
        )

    def test_assemble_megahit_parallel_paired(self):
        input_files = self.get_data_path("formatted-reads/paired-end")
        _input = SingleLanePerSamplePairedEndFastqDirFmt(input_files, mode="r")
        samples = Artifact.import_data(
            "SampleData[PairedEndSequencesWithQuality]", _input
        )

        with self.test_config:
            (out,) = self.assemble_megahit.parallel(samples)._result()

        out.validate()
        self.assertIs(out.format, ContigSequencesDirFmt)

    def test_assemble_megahit_parallel_single(self):
        input_files = self.get_data_path("formatted-reads/single-end")
        _input = SingleLanePerSampleSingleEndFastqDirFmt(input_files, mode="r")
        samples = Artifact.import_data("SampleData[SequencesWithQuality]", _input)

        with self.test_config:
            (out,) = self.assemble_megahit.parallel(samples)._result()

        out.validate()
        self.assertIs(out.format, ContigSequencesDirFmt)

    @parameterized.expand([("shortuuid",), ("uuid3",), ("uuid4",), ("uuid5",)])
    @patch("q2_assembly.megahit.megahit.modify_contig_ids")
    @patch("q2_assembly.megahit.megahit._process_sample")
    def test_assemble_megahit_different_uuids(self, uuid_type, p1, p2):
        input_files = self.get_data_path("reads/single-end")
        input = SingleLanePerSampleSingleEndFastqDirFmt(input_files, mode="r")

        obs = assemble_megahit_helper(
            reads=input,
            coassemble=False,
            uuid_type=uuid_type,
            common_args=self.test_params_list,
        )

        p2.assert_has_calls(
            [call(os.path.join(str(obs), f"sample1_contigs.fa"), "sample1", uuid_type), call(os.path.join(str(obs), f"sample2_contigs.fa"), "sample2", uuid_type)]
        )


if __name__ == "__main__":
    unittest.main()
