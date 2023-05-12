# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import contextlib
import json
import os
import tempfile
import unittest
from subprocess import CalledProcessError
from unittest.mock import ANY, patch

from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt
from qiime2.plugin.testing import TestPluginBase

from ..quast import _evaluate_contigs, _process_quast_arg, evaluate_contigs


class MockTempDir(tempfile.TemporaryDirectory):
    pass


class TestQuast(TestPluginBase):
    package = "q2_assembly.tests"

    def setUp(self):
        super().setUp()
        with contextlib.ExitStack() as stack:
            self._tmp = stack.enter_context(tempfile.TemporaryDirectory())
            self.addCleanup(stack.pop_all().close)

    def test_process_quast_arg_simple1(self):
        obs = _process_quast_arg("not_k_list", 123)
        exp = ["--not-k-list", "123"]
        self.assertListEqual(obs, exp)

    def test_process_quast_arg_simple2(self):
        obs = _process_quast_arg("k_list", [1, 2, 3])
        exp = ["--k-list", "1,2,3"]
        self.assertListEqual(obs, exp)

    @patch("platform.system", return_value="Darwin")
    def test_process_quast_arg_threads_many_Darwin(self, p1):
        with self.assertRaisesRegex(ValueError, "only supported on Linux"):
            _process_quast_arg("threads", 6)

    @patch("platform.system", return_value="Linux")
    def test_process_quast_arg_threads_many_Linux(self, p1):
        obs = _process_quast_arg("threads", 6)
        exp = ["--threads", "6"]
        self.assertListEqual(obs, exp)

    @patch("platform.system", return_value="")
    def test_process_quast_arg_threads_many_unknownOS(self, p1):
        with self.assertRaisesRegex(ValueError, "only supported on Linux"):
            _process_quast_arg("threads", 6)

    def test_process_quast_arg_threads_correct(self):
        obs = _process_quast_arg("threads", 1)
        exp = ["--threads", "1"]
        self.assertListEqual(obs, exp)

    def test_process_quast_arg_bool(self):
        obs = _process_quast_arg("k_bool", True)
        exp = ["--k-bool"]
        self.assertListEqual(obs, exp)

    @patch("subprocess.run")
    def test_evaluate_contigs_minimal(self, p):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        obs_samples = _evaluate_contigs(
            results_dir="some/dir",
            contigs=contigs,
            reads={},
            paired=False,
            references=None,
            common_args=["-t", "1"],
        )

        exp_command = [
            "metaquast.py",
            "-o",
            "some/dir",
            "-t",
            "1",
            os.path.join(str(contigs), "sample1_contigs.fa"),
            os.path.join(str(contigs), "sample2_contigs.fa"),
        ]
        self.assertListEqual(obs_samples, ["sample1", "sample2"])
        p.assert_called_once_with(exp_command, check=True)

    @patch("subprocess.run")
    def test_evaluate_contigs_more_params(self, p):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        obs_samples = _evaluate_contigs(
            results_dir="some/dir",
            contigs=contigs,
            reads={},
            paired=False,
            references=None,
            common_args=["-m", "10", "-t", "1"],
        )

        exp_command = [
            "metaquast.py",
            "-o",
            "some/dir",
            "-m",
            "10",
            "-t",
            "1",
            os.path.join(str(contigs), "sample1_contigs.fa"),
            os.path.join(str(contigs), "sample2_contigs.fa"),
        ]
        self.assertListEqual(obs_samples, ["sample1", "sample2"])
        p.assert_called_once_with(exp_command, check=True)

    @patch("subprocess.run")
    def test_evaluate_contigs_single_end(self, p):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        reads = {
            "sample1": {"fwd": "path/to/s1fwd", "rev": None},
            "sample2": {"fwd": "path/to/s2fwd", "rev": None},
        }
        obs_samples = _evaluate_contigs(
            results_dir="some/dir",
            contigs=contigs,
            reads=reads,
            paired=False,
            references=None,
            common_args=["-m", "10", "-t", "1"],
        )

        exp_command = [
            "metaquast.py",
            "-o",
            "some/dir",
            "-m",
            "10",
            "-t",
            "1",
            os.path.join(str(contigs), "sample1_contigs.fa"),
            os.path.join(str(contigs), "sample2_contigs.fa"),
            "--single",
            "path/to/s1fwd",
            "--single",
            "path/to/s2fwd",
        ]
        self.assertListEqual(obs_samples, ["sample1", "sample2"])
        p.assert_called_once_with(exp_command, check=True)

    @patch("subprocess.run")
    def test_evaluate_contigs_paired_end(self, p):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        reads = {
            "sample1": {"fwd": "path/to/s1fwd", "rev": "path/to/s1rev"},
            "sample2": {"fwd": "path/to/s2fwd", "rev": "path/to/s2rev"},
        }
        obs_samples = _evaluate_contigs(
            results_dir="some/dir",
            contigs=contigs,
            reads=reads,
            paired=True,
            references=None,
            common_args=["-m", "10", "-t", "1"],
        )

        exp_command = [
            "metaquast.py",
            "-o",
            "some/dir",
            "-m",
            "10",
            "-t",
            "1",
            os.path.join(str(contigs), "sample1_contigs.fa"),
            os.path.join(str(contigs), "sample2_contigs.fa"),
            "--pe1",
            "path/to/s1fwd",
            "--pe2",
            "path/to/s1rev",
            "--pe1",
            "path/to/s2fwd",
            "--pe2",
            "path/to/s2rev",
        ]
        self.assertListEqual(obs_samples, ["sample1", "sample2"])
        p.assert_called_once_with(exp_command, check=True)

    def test_evaluate_contigs_missing_rev_reads(self):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        reads = {
            "sample1": {"fwd": "path/to/s1fwd", "rev": "path/to/s1rev"},
        }
        with self.assertRaisesRegex(
            Exception, r".*reverse reads \(1\) does not match.*contig files \(2\).*"
        ):
            _ = _evaluate_contigs(
                results_dir="some/dir",
                contigs=contigs,
                reads=reads,
                paired=True,
                references=None,
                common_args=["-m", "10", "-t", "1"],
            )

    def test_evaluate_contigs_non_matching_samples(self):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        reads = {
            "sample1": {"fwd": "path/to/s1fwd", "rev": "path/to/s1rev"},
            "sample3": {"fwd": "path/to/s3fwd", "rev": "path/to/s3rev"},
        }
        with self.assertRaisesRegex(
            Exception, "Some samples are missing from the reads file."
        ):
            _ = _evaluate_contigs(
                results_dir="some/dir",
                contigs=contigs,
                reads=reads,
                paired=True,
                references=None,
                common_args=["-m", "10", "-t", "1"],
            )

    @patch(
        "subprocess.run", side_effect=CalledProcessError(returncode=123, cmd="some cmd")
    )
    def test_evaluate_contigs_with_error(self, p):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        reads = {
            "sample1": {"fwd": "path/to/s1fwd", "rev": "path/to/s1rev"},
            "sample2": {"fwd": "path/to/s2fwd", "rev": "path/to/s2rev"},
        }
        with self.assertRaisesRegex(
            Exception, r"An error.*while running QUAST.*code 123"
        ):
            _ = _evaluate_contigs(
                results_dir="some/dir",
                contigs=contigs,
                reads=reads,
                paired=True,
                references=None,
                common_args=["-m", "10", "-t", "1"],
            )

    @patch("platform.system", return_value="Linux")
    @patch("q2_assembly.quast._evaluate_contigs", return_value=["sample1", "sample2"])
    @patch("q2_assembly.quast._fix_html_reports", return_value=None)
    @patch("q2templates.render")
    @patch("tempfile.TemporaryDirectory")
    def test_evaluate_contigs_action_no_reads(self, p1, p2, p3, p4, p5):
        test_temp_dir = MockTempDir()
        os.mkdir(os.path.join(test_temp_dir.name, "results"))
        p1.return_value = test_temp_dir
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")

        evaluate_contigs(
            output_dir=self._tmp,
            contigs=contigs,
            reads=None,
            min_contig=150,
            threads=5,
            contig_thresholds=[10, 20],
        )

        p4.assert_called_once_with(
            os.path.join(test_temp_dir.name, "results"),
            contigs,
            {},
            False,
            None,
            ["--min-contig", "150", "--threads", "5", "--contig-thresholds", "10,20"],
        )
        p3.assert_called_once_with(os.path.join(test_temp_dir.name, "results"))

        exp_context = {
            "tabs": [
                {"title": "QC report", "url": "index.html"},
                {"title": "Contig browser", "url": "q2_icarus.html"},
            ],
            "samples": json.dumps(["sample1", "sample2"]),
        }
        p2.assert_called_once_with(ANY, self._tmp, context=exp_context)

    @patch("q2_assembly.quast._evaluate_contigs", return_value=["sample1", "sample2"])
    @patch("q2_assembly.quast._fix_html_reports", return_value=None)
    @patch("q2templates.render")
    @patch("tempfile.TemporaryDirectory")
    def test_evaluate_contigs_action_single_end(self, p1, p2, p3, p4):
        test_temp_dir = MockTempDir()
        os.mkdir(os.path.join(test_temp_dir.name, "results"))
        p1.return_value = test_temp_dir
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        reads = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path("reads/single-end"), "r"
        )

        evaluate_contigs(
            output_dir=self._tmp,
            contigs=contigs,
            reads=reads,
            min_contig=150,
            threads=1,
        )

        exp_reads_dict = {
            "sample1": {
                "fwd": self.get_data_path("reads/single-end/reads1_R1.fastq.gz"),
                "rev": None,
            },
            "sample2": {
                "fwd": self.get_data_path("reads/single-end/reads2_R1.fastq.gz"),
                "rev": None,
            },
        }
        p4.assert_called_once_with(
            os.path.join(test_temp_dir.name, "results"),
            contigs,
            exp_reads_dict,
            False,
            None,
            ["--min-contig", "150", "--threads", "1"],
        )
        p3.assert_called_once_with(os.path.join(test_temp_dir.name, "results"))

        exp_context = {
            "tabs": [
                {"title": "QC report", "url": "index.html"},
                {"title": "Contig browser", "url": "q2_icarus.html"},
            ],
            "samples": json.dumps(["sample1", "sample2"]),
        }
        p2.assert_called_once_with(ANY, self._tmp, context=exp_context)

    @patch("q2_assembly.quast._evaluate_contigs", return_value=["sample1", "sample2"])
    @patch("q2_assembly.quast._fix_html_reports", return_value=None)
    @patch("q2templates.render")
    @patch("tempfile.TemporaryDirectory")
    def test_evaluate_contigs_action_paired_end(self, p1, p2, p3, p4):
        test_temp_dir = MockTempDir()
        os.mkdir(os.path.join(test_temp_dir.name, "results"))
        p1.return_value = test_temp_dir
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        reads = SingleLanePerSamplePairedEndFastqDirFmt(
            self.get_data_path("reads/paired-end"), "r"
        )

        evaluate_contigs(
            output_dir=self._tmp,
            contigs=contigs,
            reads=reads,
            min_contig=150,
            threads=1,
        )

        exp_reads_dict = {
            "sample1": {
                "fwd": self.get_data_path("reads/paired-end/reads1_R1.fastq.gz"),
                "rev": self.get_data_path("reads/paired-end/reads1_R2.fastq.gz"),
            },
            "sample2": {
                "fwd": self.get_data_path("reads/paired-end/reads2_R1.fastq.gz"),
                "rev": self.get_data_path("reads/paired-end/reads2_R2.fastq.gz"),
            },
        }
        p4.assert_called_once_with(
            os.path.join(test_temp_dir.name, "results"),
            contigs,
            exp_reads_dict,
            True,
            None,
            ["--min-contig", "150", "--threads", "1"],
        )
        p3.assert_called_once_with(os.path.join(test_temp_dir.name, "results"))

        exp_context = {
            "tabs": [
                {"title": "QC report", "url": "index.html"},
                {"title": "Contig browser", "url": "q2_icarus.html"},
            ],
            "samples": json.dumps(["sample1", "sample2"]),
        }
        p2.assert_called_once_with(ANY, self._tmp, context=exp_context)


if __name__ == "__main__":
    unittest.main()
