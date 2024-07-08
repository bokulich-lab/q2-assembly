# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
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

    def test_parse_columns_some_cols(self):
        quast_results_path = self.get_data_path("quast-results")
        transposed_report_path = os.path.join(
            quast_results_path, "transposed_report_some_cols.tsv"
        )
        transposed_report = pd.read_csv(transposed_report_path, sep="\t", header=0)
        refined_report = _parse_columns(
            transposed_report, self.some_cols_contig_thresholds
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

    def test_parse_columns_all_cols(self):
        quast_results_path = self.get_data_path("quast-results")
        transposed_report_path = os.path.join(
            quast_results_path, "transposed_report_all_cols.tsv"
        )
        transposed_report = pd.read_csv(transposed_report_path, sep="\t", header=0)
        refined_report = _parse_columns(
            transposed_report, self.all_cols_contig_thresholds
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
            transposed_report, self.all_cols_contig_thresholds
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
            transposed_report, self.all_cols_contig_thresholds
        )
        id_col = refined_report.index

        assert not any(
            "_contig" in value for value in id_col
        ), "Values in the id column contain '_contig'"


if __name__ == "__main__":
    unittest.main()
