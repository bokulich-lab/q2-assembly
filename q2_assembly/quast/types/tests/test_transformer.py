# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
import qiime2
from qiime2.plugin.testing import TestPluginBase

from q2_assembly.quast.types import QUASTResultsFormat
from q2_assembly.quast.types._transformer import _read_dataframe


class TestQUASTTransformers(TestPluginBase):
    package = "q2_assembly.quast.types.tests"

    def setUp(self):
        super().setUp()
        self.fp = self.get_data_path("quast_results.tsv")

    def test_read_dataframe(self):
        obs = _read_dataframe(self.fp)

        self.assertIsInstance(obs, pd.DataFrame)
        self.assertEqual(obs.shape, (2, 42))
        self.assertEqual(obs.index.name, "id")

    def test_result_to_dataframe_transformer(self):
        transformer = self.get_transformer(QUASTResultsFormat, pd.DataFrame)
        obs = transformer(QUASTResultsFormat(self.fp, mode="r"))
        exp = pd.read_csv(self.fp, sep="\t", header=0, index_col=0, dtype="str")
        exp.index.name = "id"

        pd.testing.assert_frame_equal(obs, exp)

    def test_dataframe_to_result_transformer(self):
        transformer = self.get_transformer(pd.DataFrame, QUASTResultsFormat)
        df = pd.read_csv(self.fp, sep="\t", header=0, index_col=False, dtype="str")
        df.index.name = "id"
        obs = transformer(df)

        obs.validate()
        self.assertIsInstance(obs, QUASTResultsFormat)

    def test_result_to_metadata_transformer(self):
        transformer = self.get_transformer(QUASTResultsFormat, qiime2.Metadata)
        obs = transformer(QUASTResultsFormat(self.fp, mode="r"))

        df = pd.read_csv(self.fp, sep="\t", header=0, index_col=0, dtype="str")
        exp = qiime2.Metadata(df)

        self.assertEqual(obs, exp)
