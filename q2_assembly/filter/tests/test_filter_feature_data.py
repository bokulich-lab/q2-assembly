# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os

import pandas as pd
import qiime2 as q2
from q2_types.per_sample_sequences import ContigSequencesDirFmt

from q2_assembly.filter.filter import _filter_contigs
from q2_assembly.filter.tests._utils import FilterContigsTestBase

# contig IDs in tests/data/pooled/pooled_contigs.fa, in file order
POOLED_SAMPLE1_IDS = [
    f"sample1:dPj72SpWBxMXAVyPtYqJG{i}" for i in ["A", 2, 3, 4, 5, 6, 7, 8, 9, "B"]
]
POOLED_SAMPLE2_IDS = [f"sample2:2Ta5AWdPGsJiCUZJ4NeP9{i}" for i in ["A", 2, 3, 4]]


class TestFilterFeatureDataContigsByContig(FilterContigsTestBase):

    def setUp(self):
        super().setUp()
        self.contigs_pooled = ContigSequencesDirFmt(self.get_data_path("pooled"), "r")
        self.contigs_coassembled = ContigSequencesDirFmt(
            self.get_data_path("coassembled"), "r"
        )

    def test_filter_pooled_whole_and_short_ids(self):
        obs = _filter_contigs(
            contigs=self.contigs_pooled,
            on="contig",
            ids=["sample1:dPj72SpWBxMXAVyPtYqJGA", "2Ta5AWdPGsJiCUZJ4NeP92"],
        )
        self.assertDictEqual(
            obs.sample_dict(), {"pooled": os.path.join(obs.path, "pooled_contigs.fa")}
        )
        self.assertListEqual(
            self._read_all_contig_ids(obs),
            ["sample1:dPj72SpWBxMXAVyPtYqJGA", "sample2:2Ta5AWdPGsJiCUZJ4NeP92"],
        )

    def test_filter_pooled_short_ids_exclude(self):
        obs = _filter_contigs(
            contigs=self.contigs_pooled,
            on="contig",
            ids=[
                "2Ta5AWdPGsJiCUZJ4NeP9A",
                "2Ta5AWdPGsJiCUZJ4NeP92",
                "2Ta5AWdPGsJiCUZJ4NeP93",
                "2Ta5AWdPGsJiCUZJ4NeP94",
            ],
            exclude_ids=True,
        )
        self.assertListEqual(self._read_all_contig_ids(obs), POOLED_SAMPLE1_IDS)

    def test_filter_coassembled_ids(self):
        obs = _filter_contigs(
            contigs=self.contigs_coassembled,
            on="contig",
            ids=["k141_0", "k141_2"],
        )
        self.assertDictEqual(
            obs.sample_dict(),
            {"coassembled": os.path.join(obs.path, "coassembled_contigs.fa")},
        )
        self.assertListEqual(self._read_all_contig_ids(obs), ["k141_2", "k141_0"])

    def test_filter_coassembled_everything(self):
        with self.assertWarnsRegex(UserWarning, "No contigs remain after filtering"):
            obs = _filter_contigs(
                contigs=self.contigs_coassembled,
                on="contig",
                ids=["k141_0"],
                exclude_ids=True,
                length_threshold=500,
            )
        self.assertDictEqual(
            obs.sample_dict(),
            {"coassembled": os.path.join(obs.path, "coassembled_contigs.fa")},
        )
        self.assertListEqual(self._read_all_contig_ids(obs), [])


class TestFilterPooledContigsBySample(FilterContigsTestBase):

    def setUp(self):
        super().setUp()
        self.pooled_contigs = ContigSequencesDirFmt(self.get_data_path("pooled"), "r")

    def test_filter_ids(self):
        obs = _filter_contigs(contigs=self.pooled_contigs, on="sample", ids=["sample1"])
        self.assertListEqual(self._read_all_contig_ids(obs), POOLED_SAMPLE1_IDS)

    def test_filter_ids_exclude(self):
        obs = _filter_contigs(
            contigs=self.pooled_contigs,
            on="sample",
            ids=["sample1"],
            exclude_ids=True,
        )
        self.assertListEqual(self._read_all_contig_ids(obs), POOLED_SAMPLE2_IDS)

    def test_filter_metadata_and_length(self):
        metadata = q2.Metadata(
            pd.DataFrame(
                data={"col1": ["yes", "no"]},
                index=pd.Index(["sample1", "sample2"], name="id"),
            )
        )
        obs = _filter_contigs(
            contigs=self.pooled_contigs,
            on="sample",
            metadata=metadata,
            where="col1='yes'",
            length_threshold=400,
        )
        self.assertListEqual(
            self._read_all_contig_ids(obs),
            ["sample1:dPj72SpWBxMXAVyPtYqJG7", "sample1:dPj72SpWBxMXAVyPtYqJG9"],
        )

    def test_filter_missing_ids(self):
        with self.assertRaisesRegex(
            ValueError, "not present in the contig data: sample3"
        ):
            _filter_contigs(contigs=self.pooled_contigs, on="sample", ids=["sample3"])

    def test_filter_metadata_with_contig_ids(self):
        metadata = q2.Metadata(
            pd.DataFrame(
                data={"col1": ["yes"]},
                index=pd.Index(["sample1:dPj72SpWBxMXAVyPtYqJGA"], name="id"),
            )
        )
        with self.assertWarnsRegex(UserWarning, "No contigs remain after filtering"):
            obs = _filter_contigs(
                contigs=self.pooled_contigs, on="sample", metadata=metadata
            )
        self.assertListEqual(self._read_all_contig_ids(obs), [])

    def test_empty_metadata_query_include(self):
        metadata = q2.Metadata(
            pd.DataFrame(
                {"col1": ["yes"]},
                index=pd.Index(["sample1"], name="id"),
            )
        )
        with self.assertWarnsRegex(UserWarning, "No contigs remain after filtering"):
            obs = _filter_contigs(
                contigs=self.pooled_contigs,
                on="sample",
                metadata=metadata,
                where="col1='no'",
            )
        self.assertDictEqual(
            obs.sample_dict(relative=True),
            self.pooled_contigs.sample_dict(relative=True),
        )
        self.assertListEqual(self._read_all_contig_ids(obs), [])

    def test_empty_metadata_query_exclude(self):
        metadata = q2.Metadata(
            pd.DataFrame(
                {"col1": ["yes"]},
                index=pd.Index(["sample1"], name="id"),
            )
        )
        obs = _filter_contigs(
            contigs=self.pooled_contigs,
            on="sample",
            metadata=metadata,
            where="col1='no'",
            exclude_ids=True,
        )
        self.assertListEqual(
            self._read_all_contig_ids(obs), POOLED_SAMPLE1_IDS + POOLED_SAMPLE2_IDS
        )


class TestFilterFeatureDataContigsPipeline(FilterContigsTestBase):

    def setUp(self):
        super().setUp()
        self.filter_contigs = self.plugin.pipelines["filter_contigs"]
        self.coassembled_contigs = q2.Artifact.import_data(
            "FeatureData[Contig % Properties('coassembly')]",
            self.get_data_path("coassembled"),
        )
        self.pooled_contigs = q2.Artifact.import_data(
            "FeatureData[Contig % Properties('pooled')]",
            self.get_data_path("pooled"),
        )
        self.sample_metadata = q2.Metadata(
            pd.DataFrame(
                {"selected": ["yes"]},
                index=pd.Index(["sample1"], name="id"),
            )
        )

    def test_pooled_contigs_on_sample(self):
        (obs,) = self.filter_contigs(
            contigs=self.pooled_contigs,
            metadata=self.sample_metadata,
            where="selected='yes'",
        )
        self.assertEqual(obs.type, self.pooled_contigs.type)
        self.assertDictEqual(
            obs.view(ContigSequencesDirFmt).sample_dict(relative=True),
            self.pooled_contigs.view(ContigSequencesDirFmt).sample_dict(relative=True),
        )
        self.assertListEqual(self._read_all_contig_ids(obs), POOLED_SAMPLE1_IDS)

    def test_coassembled_contigs_on_contig(self):
        (obs,) = self.filter_contigs(
            contigs=self.coassembled_contigs, on="contig", ids=["k141_0"]
        )

        self.assertEqual(obs.type, self.coassembled_contigs.type)
        self.assertListEqual(self._read_all_contig_ids(obs), ["k141_0"])

    def test_pooled_contigs_on_contig(self):
        (obs,) = self.filter_contigs(
            contigs=self.pooled_contigs,
            on="contig",
            ids=["dPj72SpWBxMXAVyPtYqJG2"],
        )

        self.assertEqual(obs.type, self.pooled_contigs.type)
        self.assertListEqual(
            self._read_all_contig_ids(obs), ["sample1:dPj72SpWBxMXAVyPtYqJG2"]
        )

    def test_coassembled_contigs_on_sample_not_allowed(self):
        # on="sample" is the default
        with self.assertRaisesRegex(ValueError, "No solution for inputs"):
            self.filter_contigs(contigs=self.coassembled_contigs, ids=["k141_0"])

    def test_pooled_contigs_remove_empty_not_allowed(self):
        with self.assertRaisesRegex(ValueError, "No solution for inputs"):
            self.filter_contigs(
                contigs=self.pooled_contigs, ids=["sample1"], remove_empty=True
            )

    def test_coassembled_contigs_remove_empty_not_allowed(self):
        with self.assertRaisesRegex(ValueError, "No solution for inputs"):
            self.filter_contigs(
                contigs=self.coassembled_contigs,
                on="contig",
                ids=["k141_0"],
                remove_empty=True,
            )
