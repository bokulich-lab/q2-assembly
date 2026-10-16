# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import contextlib
import os
import tempfile
import unittest

import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_assembly.quast.report import initialize_optional_cols_map
from q2_assembly.quast.types import QUASTResultsFormat

from ..quast import _parse_columns


class MockTempDir(tempfile.TemporaryDirectory):
    pass


class TestQuastUtils(TestPluginBase):
    package = "q2_assembly.tests"

    def setUp(self):
        super().setUp()
        with contextlib.ExitStack() as stack:
            self._tmp = stack.enter_context(tempfile.TemporaryDirectory())
            self.addCleanup(stack.pop_all().close)
        self.all_cols_contig_thresholds = [1000, 5000, 10000, 25000, 50000]
        self.some_cols_contig_thresholds = [1000, 5000]
        self.sample_ids = {"sample1_contigs": "sample1", "sample2_contigs": "sample2"}

    def test_parse_columns_some_cols(self):
        quast_results_path = self.get_data_path("quast-results")
        transposed_report_path = os.path.join(
            quast_results_path, "transposed_report_some_cols.tsv"
        )
        transposed_report = pd.read_csv(transposed_report_path, sep="\t", header=0)
        refined_report = _parse_columns(
            transposed_report, self.some_cols_contig_thresholds, self.sample_ids
        )
        refined_reports_cols = refined_report.columns.tolist()
        true_columns = QUASTResultsFormat.HEADER + [
            "total_length_1000",
            "no_contigs_1000",
            "no_contigs_5000",
            "reference_length",
            "total_length_5000",
            "total_length_0",
        ]

        true_columns.remove("id")

        assert set(true_columns) == set(refined_reports_cols)

    def test_parse_columns_preserves_sample_ids(self):
        report = pd.read_csv(
            self.get_data_path("quast-results/transposed_report_sample_ids.tsv"),
            sep="\t",
        )
        result = _parse_columns(
            report,
            [],
            {
                "patient_contigs_A": "patient_contigs_A",
                "patient_A_contigs": "patient_A",
            },
        )
        self.assertListEqual(result.index.tolist(), ["patient_contigs_A", "patient_A"])

    def test_parse_columns_unknown_label(self):
        report = pd.read_csv(
            self.get_data_path("quast-results/transposed_report_sample_ids.tsv"),
            sep="\t",
        )
        with self.assertRaises(KeyError):
            _parse_columns(report, [], {})

    def test_parse_columns_all_cols(self):
        quast_results_path = self.get_data_path("quast-results")
        transposed_report_path = os.path.join(
            quast_results_path, "transposed_report_all_cols.tsv"
        )
        transposed_report = pd.read_csv(transposed_report_path, sep="\t", header=0)
        refined_report = _parse_columns(
            transposed_report, self.all_cols_contig_thresholds, self.sample_ids
        )
        refined_reports_cols = refined_report.columns.tolist()
        optional_cols_map = initialize_optional_cols_map(
            self.all_cols_contig_thresholds
        )
        true_columns = QUASTResultsFormat.HEADER + list(optional_cols_map.values())
        true_columns.remove("id")

        assert set(refined_reports_cols) == set(true_columns)

    def test_parse_columns_mandatory_cols(self):
        quast_results_path = self.get_data_path("quast-results")
        transposed_report_path = os.path.join(
            quast_results_path, "transposed_report_mandatory_cols.tsv"
        )
        transposed_report = pd.read_csv(transposed_report_path, sep="\t", header=0)
        refined_report = _parse_columns(
            transposed_report, self.all_cols_contig_thresholds, self.sample_ids
        )
        refined_reports_cols = refined_report.columns.tolist()
        true_columns = QUASTResultsFormat.HEADER.copy()
        true_columns.remove("id")

        assert set(refined_reports_cols) == set(true_columns)

    def test_parse_columns_check_id_col(self):
        quast_results_path = self.get_data_path("quast-results")
        transposed_report_path = os.path.join(
            quast_results_path, "transposed_report_all_cols.tsv"
        )
        transposed_report = pd.read_csv(transposed_report_path, sep="\t", header=0)
        refined_report = _parse_columns(
            transposed_report, self.all_cols_contig_thresholds, self.sample_ids
        )
        self.assertListEqual(refined_report.index.tolist(), ["sample1", "sample2"])


if __name__ == "__main__":
    unittest.main()
