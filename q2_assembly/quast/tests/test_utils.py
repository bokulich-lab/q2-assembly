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

from q2_assembly.quast.quast_report.report_headers import (
    OPTIONAL_COLS_MAP,
    RESHUFFLED_COLUMNS,
)

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

    def test_parse_columns_some_cols(self):
        quast_results_path = self.get_data_path("quast-results")
        transposed_report_path = os.path.join(
            quast_results_path, "transposed_report_some_cols.tsv"
        )
        transposed_report = pd.read_csv(transposed_report_path, sep="\t", header=0)
        refined_report = _parse_columns(transposed_report)
        refined_reports_cols = refined_report.columns.tolist()
        true_columns = RESHUFFLED_COLUMNS + [
            "total_length_1000",
            "no_contigs_1000",
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
        refined_report = _parse_columns(transposed_report)
        refined_reports_cols = refined_report.columns.tolist()
        true_columns = RESHUFFLED_COLUMNS + list(OPTIONAL_COLS_MAP.values())
        true_columns.remove("id")

        assert set(refined_reports_cols) == set(true_columns)

    def test_parse_columns_mandatory_cols(self):
        quast_results_path = self.get_data_path("quast-results")
        transposed_report_path = os.path.join(
            quast_results_path, "transposed_report_mandatory_cols.tsv"
        )
        transposed_report = pd.read_csv(transposed_report_path, sep="\t", header=0)
        refined_report = _parse_columns(transposed_report)
        refined_reports_cols = refined_report.columns.tolist()
        true_columns = RESHUFFLED_COLUMNS.copy()
        true_columns.remove("id")

        assert set(refined_reports_cols) == set(true_columns)

    def test_parse_columns_check_id_col(self):
        quast_results_path = self.get_data_path("quast-results")
        transposed_report_path = os.path.join(
            quast_results_path, "transposed_report_all_cols.tsv"
        )
        transposed_report = pd.read_csv(transposed_report_path, sep="\t", header=0)
        refined_report = _parse_columns(transposed_report)
        id_col = refined_report.index

        assert not any(
            "_contig" in value for value in id_col
        ), "Values in the id column contain '_contig'"


if __name__ == "__main__":
    unittest.main()
