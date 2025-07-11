# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import contextlib
import json
import os
import shutil
import tempfile
import unittest
import warnings
from filecmp import dircmp
from subprocess import CalledProcessError
from unittest.mock import ANY, MagicMock, Mock, call, mock_open, patch
from zipfile import ZipFile

import pandas as pd
from q2_types.feature_data import DNAFASTAFormat
from q2_types.genome_data import GenomeSequencesDirectoryFormat
from q2_types.per_sample_sequences import (
    BAMDirFmt,
    ContigSequencesDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from qiime2 import Artifact
from qiime2.core.exceptions import ValidationError
from qiime2.plugin.testing import TestPluginBase
from qiime2.sdk import Context

from ..quast import (
    _create_tabular_results,
    _evaluate_quast,
    _move_references,
    _process_quast_arg,
    _split_reference,
    _visualize_quast,
    _zip_additional_reports,
    _zip_dir,
    evaluate_quast,
)


class MockTempDir(tempfile.TemporaryDirectory):
    pass


class TestQuast(TestPluginBase):
    package = "q2_assembly.tests"

    def setUp(self):
        super().setUp()
        with contextlib.ExitStack() as stack:
            self._tmp = stack.enter_context(tempfile.TemporaryDirectory())
            self.addCleanup(stack.pop_all().close)

    def assert_is_dataframe(self, shape=None):
        class Matcher:
            def __eq__(self, other):
                is_df = isinstance(other, pd.DataFrame)
                if shape is not None:
                    if not other.shape == shape:
                        print(f"Shapes do not match: {other.shape} != {shape}")
                    return is_df and other.shape == shape
                else:
                    return is_df

        return Matcher()

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

    def test_split_reference(self):
        ref = DNAFASTAFormat(self.get_data_path("references/ref1.fasta"), "r")
        refs_dir = "some/dir"
        m = mock_open()

        with patch("builtins.open", m):
            obs = _split_reference(ref, refs_dir)
        exp1, exp2 = ["some/dir/ref1.1.fasta", "some/dir/ref1.2.fasta"]

        self.assertListEqual(obs, [exp1, exp2])
        self.assertListEqual(m.call_args_list, [call(exp1, "w"), call(exp2, "w")])
        m().write.assert_has_calls(
            [call(">ref1.1\nATGCATCG"), call(">ref1.2\nGTCATCAT")]
        )

    @patch("subprocess.run")
    def test_evaluate_quast_minimal(self, p):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        obs_samples = _evaluate_quast(
            results_dir="some/dir",
            contigs=contigs,
            reads={},
            paired=False,
            references=None,
            alignment_maps=None,
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
    def test_evaluate_quast_more_params(self, p):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        obs_samples = _evaluate_quast(
            results_dir="some/dir",
            contigs=contigs,
            reads={},
            paired=False,
            references=None,
            alignment_maps=None,
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
    def test_evaluate_quast_more_params_memory_efficient(self, p):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        obs_samples = _evaluate_quast(
            results_dir="some/dir",
            contigs=contigs,
            reads={},
            paired=False,
            references=None,
            alignment_maps=None,
            common_args=["-m", "10", "-t", "1", "--memory-efficient"],
        )

        exp_command = [
            "metaquast.py",
            "-o",
            "some/dir",
            "-m",
            "10",
            "-t",
            "1",
            "--memory-efficient",
            os.path.join(str(contigs), "sample1_contigs.fa"),
            os.path.join(str(contigs), "sample2_contigs.fa"),
        ]
        self.assertListEqual(obs_samples, ["sample1", "sample2"])
        p.assert_called_once_with(exp_command, check=True)

    @patch("subprocess.run")
    def test_evaluate_quast_with_map(self, p):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        alignment_maps = BAMDirFmt(self.get_data_path("alignment_map"), "r")
        obs_samples = _evaluate_quast(
            results_dir="some/dir",
            contigs=contigs,
            reads=None,
            paired=False,
            references=None,
            alignment_maps=alignment_maps,
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
            "--bam",
            f"{os.path.join(str(alignment_maps), 'sample1_alignment.bam')},"
            f"{os.path.join(str(alignment_maps), 'sample2_alignment.bam')}",
        ]
        self.assertListEqual(obs_samples, ["sample1", "sample2"])
        p.assert_called_once_with(exp_command, check=True)

    @patch("subprocess.run")
    def test_evaluate_quast_with_map_and_reads(self, p):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        alignment_maps = BAMDirFmt(self.get_data_path("alignment_map"), "r")
        reads = {
            "sample1": {"fwd": "path/to/s1fwd", "rev": None},
            "sample2": {"fwd": "path/to/s2fwd", "rev": None},
        }
        obs_samples = _evaluate_quast(
            results_dir="some/dir",
            contigs=contigs,
            reads=reads,
            paired=False,
            references=None,
            alignment_maps=alignment_maps,
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
            "--bam",
            f"{os.path.join(str(alignment_maps), 'sample1_alignment.bam')},"
            f"{os.path.join(str(alignment_maps), 'sample2_alignment.bam')}",
        ]

        self.assertListEqual(obs_samples, ["sample1", "sample2"])
        p.assert_called_once_with(exp_command, check=True)

    @patch("subprocess.run")
    def test_evaluate_quast_single_end(self, p):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        reads = {
            "sample1": {"fwd": "path/to/s1fwd", "rev": None},
            "sample2": {"fwd": "path/to/s2fwd", "rev": None},
        }
        obs_samples = _evaluate_quast(
            results_dir="some/dir",
            contigs=contigs,
            reads=reads,
            paired=False,
            references=None,
            alignment_maps=None,
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
    def test_evaluate_quast_paired_end(self, p):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        reads = {
            "sample1": {"fwd": "path/to/s1fwd", "rev": "path/to/s1rev"},
            "sample2": {"fwd": "path/to/s2fwd", "rev": "path/to/s2rev"},
        }
        obs_samples = _evaluate_quast(
            results_dir="some/dir",
            contigs=contigs,
            reads=reads,
            paired=True,
            references=None,
            alignment_maps=None,
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

    def test_evaluate_quast_missing_rev_reads(self):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        reads = {
            "sample1": {"fwd": "path/to/s1fwd", "rev": "path/to/s1rev"},
        }
        with self.assertRaisesRegex(
            Exception, r".*reverse reads \(1\) does not match.*contig files \(2\).*"
        ):
            _ = _evaluate_quast(
                results_dir="some/dir",
                contigs=contigs,
                reads=reads,
                paired=True,
                references=None,
                alignment_maps=None,
                common_args=["-m", "10", "-t", "1"],
            )

    def test_evaluate_quast_non_matching_samples(self):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        reads = {
            "sample1": {"fwd": "path/to/s1fwd", "rev": "path/to/s1rev"},
            "sample3": {"fwd": "path/to/s3fwd", "rev": "path/to/s3rev"},
        }
        with self.assertRaisesRegex(
            Exception, "Some samples are missing from the reads file."
        ):
            _ = _evaluate_quast(
                results_dir="some/dir",
                contigs=contigs,
                reads=reads,
                paired=True,
                references=None,
                alignment_maps=None,
                common_args=["-m", "10", "-t", "1"],
            )

    @patch(
        "subprocess.run", side_effect=CalledProcessError(returncode=123, cmd="some cmd")
    )
    def test_evaluate_quast_with_error(self, p):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        reads = {
            "sample1": {"fwd": "path/to/s1fwd", "rev": "path/to/s1rev"},
            "sample2": {"fwd": "path/to/s2fwd", "rev": "path/to/s2rev"},
        }
        with self.assertRaisesRegex(
            Exception, r"An error.*while running QUAST.*code 123"
        ):
            _ = _evaluate_quast(
                results_dir="some/dir",
                contigs=contigs,
                reads=reads,
                paired=True,
                references=None,
                alignment_maps=None,
                common_args=["-m", "10", "-t", "1"],
            )

    @patch("subprocess.run")
    def test_evaluate_quast_with_refs(self, p1):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        refs = GenomeSequencesDirectoryFormat(self.get_data_path("references"), "r")

        exp_refs = [
            os.path.join(str(refs.path), "ref1.fasta"),
            os.path.join(str(refs.path), "ref2.fasta"),
        ]

        obs_samples = _evaluate_quast(
            results_dir="some/dir",
            contigs=contigs,
            reads={},
            paired=False,
            references=refs,
            alignment_maps=None,
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
            "-r",
            exp_refs[0],
            "-r",
            exp_refs[1],
        ]
        self.assertListEqual(obs_samples, ["sample1", "sample2"])
        p1.assert_called_once_with(exp_command, check=True)

    @patch("q2_assembly.quast._create_tabular_results")
    @patch("platform.system", return_value="Linux")
    @patch("q2_assembly.quast._evaluate_quast", return_value=["sample1", "sample2"])
    @patch("q2_assembly.quast._fix_html_reports", return_value=None)
    @patch("q2templates.render")
    @patch("tempfile.TemporaryDirectory")
    def test_visualize_quast_action_no_reads(self, p1, p2, p3, p4, p5, p6):
        test_temp_dir = MockTempDir()
        os.mkdir(os.path.join(test_temp_dir.name, "results"))
        p1.return_value = test_temp_dir
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        p6.return_value = pd.DataFrame({"col1": [1, 2, 3], "col2": [4, 5, 6]})
        mock_to_csv = Mock()
        p6.return_value.to_csv = mock_to_csv

        _visualize_quast(
            output_dir=self._tmp,
            contigs=contigs,
            reads=None,
            min_contig=150,
            threads=5,
            k_mer_size=101,
            contig_thresholds=[10, 20],
            genomes_dir=test_temp_dir.name,
        )

        p4.assert_called_once_with(
            os.path.join(test_temp_dir.name, "results"),
            contigs,
            {},
            False,
            None,
            None,
            [
                "--min-contig",
                "150",
                "--threads",
                "5",
                "--k-mer-size",
                "101",
                "--contig-thresholds",
                "10,20",
                "--min-alignment",
                "65",
                "--min-identity",
                "90.0",
                "--ambiguity-usage",
                "one",
                "--ambiguity-score",
                "0.99",
            ],
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

    @patch("q2_assembly.quast._create_tabular_results")
    @patch("q2_assembly.quast._evaluate_quast", return_value=["sample1", "sample2"])
    @patch("q2_assembly.quast._fix_html_reports", return_value=None)
    @patch("q2templates.render")
    @patch("tempfile.TemporaryDirectory")
    def test_visualize_quast_action_single_end(self, p1, p2, p3, p4, p5):
        test_temp_dir = MockTempDir()
        os.mkdir(os.path.join(test_temp_dir.name, "results"))
        os.mkdir(os.path.join(test_temp_dir.name, "results", "krona_charts"))
        p1.return_value = test_temp_dir
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        reads = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path("reads/single-end"), "r"
        )
        p5.return_value = pd.DataFrame({"col1": [1, 2, 3], "col2": [4, 5, 6]})
        mock_to_csv = Mock()
        p5.return_value.to_csv = mock_to_csv

        _visualize_quast(
            output_dir=self._tmp,
            contigs=contigs,
            reads=reads,
            min_contig=150,
            threads=1,
            k_mer_size=101,
            contig_thresholds=[0, 1000, 5000, 10000, 25000, 50000],
            genomes_dir=test_temp_dir.name,
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
            None,
            [
                "--min-contig",
                "150",
                "--threads",
                "1",
                "--k-mer-size",
                "101",
                "--contig-thresholds",
                "0,1000,5000,10000,25000,50000",
                "--min-alignment",
                "65",
                "--min-identity",
                "90.0",
                "--ambiguity-usage",
                "one",
                "--ambiguity-score",
                "0.99",
            ],
        )
        p3.assert_called_once_with(os.path.join(test_temp_dir.name, "results"))

        exp_context = {
            "tabs": [
                {"title": "QC report", "url": "index.html"},
                {"title": "Contig browser", "url": "q2_icarus.html"},
                {"title": "Krona charts", "url": "q2_krona_charts.html"},
            ],
            "samples": json.dumps(["sample1", "sample2"]),
        }
        p2.assert_called_once_with(ANY, self._tmp, context=exp_context)

    @patch("q2_assembly.quast._create_tabular_results")
    @patch("q2_assembly.quast._evaluate_quast", return_value=["sample1", "sample2"])
    @patch("q2_assembly.quast._fix_html_reports", return_value=None)
    @patch("q2templates.render")
    @patch("tempfile.TemporaryDirectory")
    def test_visualize_quast_action_paired_end(self, p1, p2, p3, p4, p5):
        test_temp_dir = MockTempDir()
        os.mkdir(os.path.join(test_temp_dir.name, "results"))
        p1.return_value = test_temp_dir
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        reads = SingleLanePerSamplePairedEndFastqDirFmt(
            self.get_data_path("reads/paired-end"), "r"
        )
        p5.return_value = pd.DataFrame({"col1": [1, 2, 3], "col2": [4, 5, 6]})
        mock_to_csv = Mock()
        p5.return_value.to_csv = mock_to_csv

        _visualize_quast(
            output_dir=self._tmp,
            contigs=contigs,
            reads=reads,
            min_contig=150,
            threads=1,
            k_mer_size=101,
            contig_thresholds=[0, 1000, 5000, 10000, 25000, 50000],
            genomes_dir=test_temp_dir.name,
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
            None,
            [
                "--min-contig",
                "150",
                "--threads",
                "1",
                "--k-mer-size",
                "101",
                "--contig-thresholds",
                "0,1000,5000,10000,25000,50000",
                "--min-alignment",
                "65",
                "--min-identity",
                "90.0",
                "--ambiguity-usage",
                "one",
                "--ambiguity-score",
                "0.99",
            ],
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

    @patch("q2_assembly.quast._create_tabular_results")
    @patch("q2_assembly.quast._evaluate_quast", return_value=["sample1", "sample2"])
    @patch("q2_assembly.quast._fix_html_reports", return_value=None)
    @patch("q2templates.render")
    @patch("tempfile.TemporaryDirectory")
    def test_evaluate_quast_action_paired_end_no_icarus(self, p1, p2, p3, p4, p5):
        test_temp_dir = MockTempDir()
        os.mkdir(os.path.join(test_temp_dir.name, "results"))
        p1.return_value = test_temp_dir
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        reads = SingleLanePerSamplePairedEndFastqDirFmt(
            self.get_data_path("reads/paired-end"), "r"
        )

        _visualize_quast(
            output_dir=self._tmp,
            contigs=contigs,
            reads=reads,
            min_contig=150,
            threads=1,
            k_mer_size=101,
            contig_thresholds=[0, 1000, 5000, 10000, 25000, 50000],
            no_icarus=True,
            genomes_dir=test_temp_dir.name,
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
            None,
            [
                "--min-contig",
                "150",
                "--threads",
                "1",
                "--k-mer-size",
                "101",
                "--contig-thresholds",
                "0,1000,5000,10000,25000,50000",
                "--min-alignment",
                "65",
                "--min-identity",
                "90.0",
                "--ambiguity-usage",
                "one",
                "--ambiguity-score",
                "0.99",
                "--no-icarus",
            ],
        )
        p3.assert_called_once_with(os.path.join(test_temp_dir.name, "results"))

        exp_context = {
            "tabs": [
                {"title": "QC report", "url": "index.html"},
            ],
            "samples": json.dumps(["sample1", "sample2"]),
        }
        p2.assert_called_once_with(ANY, self._tmp, context=exp_context)

    @patch("pandas.read_csv")
    @patch("q2_assembly.quast._parse_columns")
    def test_create_tabular_results(self, p1, p2):
        report_path = os.path.join(self.temp_dir.name, "transposed_report.tsv")
        mock_df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        p2.return_value = mock_df

        _ = _create_tabular_results(self.temp_dir.name, [1000, 5000, 25000, 50000])

        p2.assert_called_once_with(report_path, sep="\t", header=0)
        p1.assert_called_once_with(mock_df, [1000, 5000, 25000, 50000])

    @patch("pandas.read_csv")
    @patch("q2_assembly.quast._parse_columns")
    def test_create_tabular_results_in_combined_subdir(self, p1, p2):
        report_path = os.path.join(
            self.temp_dir.name, "combined_reference", "transposed_report.tsv"
        )
        os.makedirs(os.path.join(self.temp_dir.name, "combined_reference"))
        mock_df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        p2.return_value = mock_df

        _ = _create_tabular_results(self.temp_dir.name, [1000, 5000, 25000, 50000])

        p2.assert_called_once_with(report_path, sep="\t", header=0)
        p1.assert_called_once_with(mock_df, [1000, 5000, 25000, 50000])

    def test_evaluate_quast_pipeline_with_refs(self):
        def _copy_references(refs_dir, tmp):
            new_refs_dir = os.path.join(tmp, "vis_files", "quast_downloaded_references")
            os.makedirs(new_refs_dir, exist_ok=True)
            for ref in os.listdir(refs_dir):
                shutil.copy(
                    src=os.path.join(refs_dir, ref), dst=os.path.join(new_refs_dir, ref)
                )

        contigs = Artifact.import_data(
            "SampleData[Contigs]", self.get_data_path("contigs")
        )

        with tempfile.TemporaryDirectory() as tmp, patch(
            "tempfile.TemporaryDirectory"
        ) as mock_temp_dir:
            mock_temp_dir.return_value.__enter__.return_value = tmp
            action, export_data, make_artifact, mock_ctx = self.create_mock_ctx()

            refs_dir = self.get_data_path("references")
            _copy_references(refs_dir, tmp)

            refs = Artifact.import_data("GenomeData[DNASequence]", refs_dir)
            tab_report, _, _ = evaluate_quast(
                ctx=mock_ctx, contigs=contigs, references=refs
            )
            make_artifact.assert_called_once_with(
                "QUASTResults", self.assert_is_dataframe((2, 43))
            )

        tab_report.validate()
        self.assertEqual(action.call_args[0][0], contigs)
        export_data.assert_called_once_with(os.path.join(tmp, "vis_files"))

    def create_mock_ctx(self):
        export_data = MagicMock(
            side_effect=lambda x: self._copy(os.path.join(x, "quast_data"))
        )
        action = MagicMock(return_value=(MagicMock(export_data=export_data),))
        ctx = Context()
        ctx._scope = MagicMock()
        make_artifact = MagicMock(side_effect=lambda *args: ctx.make_artifact(*args))
        mock_ctx = MagicMock(
            spec=Context,
            get_action=MagicMock(return_value=action),
            make_artifact=make_artifact,
        )
        return action, export_data, make_artifact, mock_ctx

    def _copy(self, dst_dir):
        os.makedirs(dst_dir, exist_ok=True)
        shutil.copy2(
            src=os.path.join(
                self.get_data_path("quast-results"), "enhanced_tabular_results.tsv"
            ),
            dst=os.path.join(dst_dir, "quast_results.tsv"),
        )

    def test_evaluate_quast_pipeline_no_refs_normal(self):
        contigs = Artifact.import_data(
            "SampleData[Contigs]", self.get_data_path("contigs")
        )

        with tempfile.TemporaryDirectory() as tmp, patch(
            "tempfile.TemporaryDirectory"
        ) as mock_temp_dir:
            mock_temp_dir.return_value.__enter__.return_value = tmp
            action, export_data, make_artifact, mock_ctx = self.create_mock_ctx()

            with patch(
                "q2_assembly.quast.GenomeSequencesDirectoryFormat"
            ) as MockGenomeSequencesDirectoryFormat:
                MockGenomeSequencesDirectoryFormat.return_value = (
                    GenomeSequencesDirectoryFormat(
                        path=self.get_data_path("references"), mode="r"
                    )
                )
                tab_report, _, _ = evaluate_quast(ctx=mock_ctx, contigs=contigs)
                make_artifact.assert_has_calls(
                    [
                        call("QUASTResults", self.assert_is_dataframe((2, 43))),
                        call("GenomeData[DNASequence]", ANY),
                    ]
                )

        tab_report.validate()
        self.assertEqual(action.call_args[0][0], contigs)
        export_data.assert_called_once_with(os.path.join(tmp, "vis_files"))

    def test_evaluate_quast_pipeline_no_refs_no_downloaded_genomes(self):
        contigs = Artifact.import_data(
            "SampleData[Contigs]", self.get_data_path("contigs")
        )

        with tempfile.TemporaryDirectory() as tmp, patch(
            "tempfile.TemporaryDirectory"
        ) as mock_temp_dir:
            mock_temp_dir.return_value.__enter__.return_value = tmp
            action, export_data, make_artifact, mock_ctx = self.create_mock_ctx()

            with patch(
                "q2_assembly.quast.GenomeSequencesDirectoryFormat"
            ) as MockGenomeSequencesDirectoryFormat:

                with warnings.catch_warnings(record=True) as w:
                    MockGenomeSequencesDirectoryFormat.return_value = (
                        GenomeSequencesDirectoryFormat()
                    )
                    tab_report, _, _ = evaluate_quast(ctx=mock_ctx, contigs=contigs)
                    make_artifact.assert_has_calls(
                        [
                            call(
                                "QUASTResults",
                                self.assert_is_dataframe((2, 43)),
                            ),
                            call("GenomeData[DNASequence]", ANY),
                        ]
                    )
                    self.assertTrue(
                        any(
                            "QUAST did not download any genomes." in str(msg.message)
                            for msg in w
                        )
                    )

            tab_report.validate()
            self.assertEqual(action.call_args[0][0], contigs)
            export_data.assert_called_once_with(os.path.join(tmp, "vis_files"))

    def test_evaluate_quast_pipeline_no_refs_corrupt_files(self):
        contigs = Artifact.import_data(
            "SampleData[Contigs]", self.get_data_path("contigs")
        )

        with tempfile.TemporaryDirectory() as tmp, patch(
            "tempfile.TemporaryDirectory"
        ) as mock_temp_dir:
            mock_temp_dir.return_value.__enter__.return_value = tmp
            action, export_data, make_artifact, mock_ctx = self.create_mock_ctx()
            with patch(
                "q2_assembly.quast.GenomeSequencesDirectoryFormat"
            ) as MockGenomeSequencesDirectoryFormat:
                with warnings.catch_warnings(record=True) as w:
                    # corrupt file
                    genomes_dir = os.path.join(tmp, "genomes_dir")
                    os.makedirs(genomes_dir, exist_ok=True)
                    corrupt_file = os.path.join(genomes_dir, "corrupt.fasta")
                    with open(corrupt_file, "w") as f:
                        f.write("seq1")
                        f.write("ACGTXCGTA")

                    mocked_dir = GenomeSequencesDirectoryFormat(
                        path=genomes_dir, mode="r"
                    )
                    MockGenomeSequencesDirectoryFormat.side_effect = [
                        mocked_dir,
                        GenomeSequencesDirectoryFormat(),
                    ]

                    tab_report, _, _ = evaluate_quast(ctx=mock_ctx, contigs=contigs)

                    self.assertRaises(ValidationError)
                    self.assertTrue(
                        any(
                            "There was a problem with "
                            "the genome files " in str(msg.message)
                            for msg in w
                        )
                    )

        make_artifact.assert_has_calls(
            [
                call("QUASTResults", self.assert_is_dataframe((2, 43))),
                call("GenomeData[DNASequence]", ANY),
            ]
        )

        tab_report.validate()
        self.assertEqual(action.call_args[0][0], contigs)
        export_data.assert_called_once_with(os.path.join(tmp, "vis_files"))

    def test_move_references(self):
        with tempfile.TemporaryDirectory() as tmp:
            refs_dir = self.get_data_path("references")
            os.makedirs(os.path.join(tmp, "references1"), exist_ok=True)
            tmp_ref_dir = os.path.join(tmp, "references1")
            empty_fp = os.path.join(tmp_ref_dir, "empty.txt")

            for file in os.listdir(refs_dir):
                shutil.copy(src=os.path.join(refs_dir, file), dst=tmp_ref_dir)

            with open(empty_fp, "w"):
                pass
            # tmp file should have the empty file
            # and the 2 reference files
            self.assertTrue(len(os.listdir(tmp_ref_dir)) == 3)
            with tempfile.TemporaryDirectory() as tmp2:
                _move_references(tmp_ref_dir, tmp2)
                files = os.listdir(tmp2)
                self.assertTrue(len(files) == 3)
                self.assertTrue(
                    ["ref1.1.fasta", "ref1.2.fasta", "ref2.1.fasta"], sorted(files)
                )

    def test_zip_dir(self):
        # Get path to test data
        data_path = self.get_data_path("zip_test_data")

        # Create tmp working directory
        with tempfile.TemporaryDirectory() as tmp:
            # Open ZipFile context and call _zip_dir
            output_filename = os.path.join(tmp, "this_is_a.zip")
            with ZipFile(output_filename, "w") as zipf:
                _zip_dir(zipf, data_path)

            with ZipFile(output_filename, "r") as zipf:
                zipf.extractall(path=os.path.join(tmp, "observed"))

            # Instantiate compare directory object
            compare = dircmp(
                a=os.path.join(data_path, "expected"), b=os.path.join(tmp, "observed")
            )

            # Should return an empty list
            self.assertFalse(compare.diff_files)

    def test_zip_additional_reports(self):
        # Get path to test data
        root_path = self.get_data_path("zip_test_data/expected")

        # Create tmp working directory
        with tempfile.TemporaryDirectory() as tmp:
            # Open ZipFile context and call _zip_dir
            output_filename = os.path.join(tmp, "this_is_a.zip")
            path_to_dirs = [os.path.join(root_path, f"folder_{i}") for i in range(3)]
            _zip_additional_reports(
                path_to_dirs=path_to_dirs, output_filename=output_filename
            )

            # Extract zip
            with ZipFile(output_filename, "r") as zipf:
                zipf.extractall(path=os.path.join(tmp, "observed"))

            # Instantiate compare directory object
            compare = dircmp(a=root_path, b=os.path.join(tmp, "observed"))

            # Should return an empty list
            self.assertFalse(compare.diff_files)


if __name__ == "__main__":
    unittest.main()
