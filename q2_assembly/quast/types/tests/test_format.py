# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.core.exceptions import ValidationError
from qiime2.plugin.testing import TestPluginBase

from q2_assembly.quast.types import QUASTResultsFormat


class TestQUASTFormats(TestPluginBase):
    package = "q2_assembly.quast.types.tests"

    def test_quast_results_format_ok(self):
        results = QUASTResultsFormat(self.get_data_path("quast_results.tsv"), mode="r")
        results.validate(level="min")
        results.validate(level="max")

    def test_quast_results_format_ok_all_cols(self):
        results = QUASTResultsFormat(
            self.get_data_path("quast_results_all_cols.tsv"), mode="r"
        )
        results.validate(level="min")
        results.validate(level="max")

    def test_quast_results_format_ok_mandatory_cols(self):
        results = QUASTResultsFormat(
            self.get_data_path("quast_results_mandatory_cols.tsv"), mode="r"
        )
        results.validate(level="min")
        results.validate(level="max")

    def test_quast_results_format_ok_some_cols(self):
        results = QUASTResultsFormat(
            self.get_data_path("quast_results_some_cols.tsv"), mode="r"
        )
        results.validate(level="min")
        results.validate(level="max")

    def test_quast_results_format_error_header(self):
        results = QUASTResultsFormat(
            self.get_data_path("quast_results_broken_header.tsv"), mode="r"
        )
        with self.assertRaisesRegex(ValidationError, "Invalid header"):
            results.validate()

    def test_quast_results_format_error_values(self):
        results = QUASTResultsFormat(
            self.get_data_path("quast_results_broken_values.tsv"), mode="r"
        )
        with self.assertRaisesRegex(ValidationError, "Line 3 has 42 columns"):
            results.validate()
