# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os

import pandas as pd
import qiime2 as q2
import skbio
from q2_types.per_sample_sequences import ContigSequencesDirFmt
from qiime2.plugin.testing import TestPluginBase

from q2_assembly.filter import filter_contigs


class TestFilterContigs(TestPluginBase):
    package = "q2_assembly.tests"

    def setUp(self):
        super().setUp()
        self.contigs = ContigSequencesDirFmt(
            self.get_data_path("contigs-with-empty"), "r"
        )
        self.metadata_df = pd.DataFrame(
            data={"col1": ["yes", "no", "yes"]},
            index=pd.Index(["sample1", "sample2", "sample3"], name="id"),
        )
        self.metadata = q2.Metadata(self.metadata_df)

    def test_filter_metadata(self):
        obs = filter_contigs(
            contigs=self.contigs, metadata=self.metadata, where="col1='yes'"
        )

        self.assertDictEqual(
            obs.sample_dict(),
            {
                "sample1": os.path.join(obs.path, "sample1_contigs.fa"),
                "sample3": os.path.join(obs.path, "sample3_contigs.fa"),
            },
        )

    def test_filter_metadata_exclude_ids(self):
        obs = filter_contigs(
            contigs=self.contigs,
            metadata=self.metadata,
            where="col1='yes'",
            exclude_ids=True,
        )

        self.assertDictEqual(
            obs.sample_dict(), {"sample2": os.path.join(obs.path, "sample2_contigs.fa")}
        )

    def test_filter_metadata_no_query_no_metadata(self):
        with self.assertRaisesRegex(
            ValueError, "At least one of the following parameters must be provided"
        ):
            filter_contigs(contigs=self.contigs)

    def test_filter_metadata_no_query(self):
        with self.assertRaisesRegex(ValueError, "A filter query must be provided"):
            filter_contigs(contigs=self.contigs, metadata=self.metadata)

    def test_filter_by_length(self):
        obs = filter_contigs(contigs=self.contigs, length_threshold=320)

        self.assertEqual(len(obs.sample_dict()), 3)

        exp_counts = (6, 1, 0)
        for (_id, fp), count in zip(obs.sample_dict().items(), exp_counts):
            with open(fp) as f:
                self.assertEqual(len(list(skbio.io.read(f, format="fasta"))), count)

    def test_filter_remove_empty(self):
        obs = filter_contigs(contigs=self.contigs, remove_empty=True)

        self.assertDictEqual(
            obs.sample_dict(),
            {
                "sample1": os.path.join(obs.path, "sample1_contigs.fa"),
                "sample2": os.path.join(obs.path, "sample2_contigs.fa"),
            },
        )

    def test_filter_by_length_and_remove_empty(self):
        obs = filter_contigs(
            contigs=self.contigs, length_threshold=400, remove_empty=True
        )

        self.assertEqual(len(obs.sample_dict()), 1)

        with open(os.path.join(obs.path, "sample1_contigs.fa")) as f:
            self.assertEqual(len(list(skbio.io.read(f, format="fasta"))), 2)

    def test_filter_everything(self):
        with self.assertRaisesRegex(ValueError, "No samples remain after filtering"):
            filter_contigs(
                contigs=self.contigs, length_threshold=1000, remove_empty=True
            )
