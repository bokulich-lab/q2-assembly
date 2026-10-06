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
import skbio
from q2_types.per_sample_sequences import ContigSequencesDirFmt
from qiime2.plugin.testing import TestPluginBase

from q2_assembly.filter import (
    _filter_contigs,
    _match_contig_ids,
    _match_sample_prefix_ids,
)

# contig IDs in tests/data/pooled/pooled_contigs.fa, in file order
POOLED_SAMPLE1_IDS = [
    f"sample1:dPj72SpWBxMXAVyPtYqJG{i}" for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 0]
]
POOLED_SAMPLE2_IDS = [f"sample2:2Ta5AWdPGsJiCUZJ4NeP9{i}" for i in [1, 2, 3, 4]]


class FilterContigsTestBase(TestPluginBase):
    package = "q2_assembly.tests"

    @staticmethod
    def _read_contig_ids(fp):
        return [c.metadata["id"] for c in skbio.io.read(fp, format="fasta")]

    def _read_all_contig_ids(self, contigs):
        if isinstance(contigs, q2.Artifact):
            contigs = contigs.view(ContigSequencesDirFmt)
        return [
            _id
            for fp in contigs.sample_dict().values()
            for _id in self._read_contig_ids(fp)
        ]


class TestFilterContigsBySample(TestPluginBase):
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

    def test_filter_sample_metadata(self):
        obs = _filter_contigs(
            contigs=self.contigs,
            on="sample",
            metadata=self.metadata,
            where="col1='yes'",
        )

        self.assertDictEqual(
            obs.sample_dict(),
            {
                "sample1": os.path.join(obs.path, "sample1_contigs.fa"),
                "sample3": os.path.join(obs.path, "sample3_contigs.fa"),
            },
        )

    def test_filter_sample_metadata_exclude_ids(self):
        obs = _filter_contigs(
            contigs=self.contigs,
            on="sample",
            metadata=self.metadata,
            where="col1='yes'",
            exclude_ids=True,
        )

        self.assertDictEqual(
            obs.sample_dict(), {"sample2": os.path.join(obs.path, "sample2_contigs.fa")}
        )

    def test_filter_sample_metadata_no_query(self):
        # When no 'where' query is provided, should filter to include all samples
        # present in the metadata (sample1, sample2, sample3)
        obs = _filter_contigs(contigs=self.contigs, on="sample", metadata=self.metadata)

        self.assertDictEqual(
            obs.sample_dict(),
            {
                "sample1": os.path.join(obs.path, "sample1_contigs.fa"),
                "sample2": os.path.join(obs.path, "sample2_contigs.fa"),
                "sample3": os.path.join(obs.path, "sample3_contigs.fa"),
            },
        )

    def test_filter_sample_metadata_subset_ids_only(self):
        # Test filtering by metadata that only contains a subset of sample IDs
        metadata_subset_df = pd.DataFrame(
            data={"col1": ["yes", "no"]},
            index=pd.Index(["sample1", "sample2"], name="id"),
        )
        metadata_subset = q2.Metadata(metadata_subset_df)

        obs = _filter_contigs(
            contigs=self.contigs, on="sample", metadata=metadata_subset
        )

        self.assertDictEqual(
            obs.sample_dict(),
            {
                "sample1": os.path.join(obs.path, "sample1_contigs.fa"),
                "sample2": os.path.join(obs.path, "sample2_contigs.fa"),
            },
        )

    def test_filter_by_length(self):
        obs = _filter_contigs(contigs=self.contigs, on="sample", length_threshold=320)

        self.assertEqual(len(obs.sample_dict()), 3)

        exp_counts = (6, 1, 0)
        for (_id, fp), count in zip(obs.sample_dict().items(), exp_counts):
            with open(fp) as f:
                self.assertEqual(len(list(skbio.io.read(f, format="fasta"))), count)

    def test_filter_remove_empty(self):
        obs = _filter_contigs(contigs=self.contigs, on="sample", remove_empty=True)

        self.assertDictEqual(
            obs.sample_dict(),
            {
                "sample1": os.path.join(obs.path, "sample1_contigs.fa"),
                "sample2": os.path.join(obs.path, "sample2_contigs.fa"),
            },
        )

    def test_filter_by_length_and_remove_empty(self):
        obs = _filter_contigs(
            contigs=self.contigs,
            on="sample",
            length_threshold=400,
            remove_empty=True,
        )

        self.assertEqual(len(obs.sample_dict()), 1)

        with open(os.path.join(obs.path, "sample1_contigs.fa")) as f:
            self.assertEqual(len(list(skbio.io.read(f, format="fasta"))), 2)

    def test_filter_everything(self):
        with self.assertRaisesRegex(ValueError, "No samples remain after filtering"):
            _filter_contigs(
                contigs=self.contigs,
                on="sample",
                length_threshold=1000,
                remove_empty=True,
            )

    def test_filter_sample_ids(self):
        obs = _filter_contigs(contigs=self.contigs, on="sample", ids=["sample2"])

        self.assertDictEqual(
            obs.sample_dict(), {"sample2": os.path.join(obs.path, "sample2_contigs.fa")}
        )

    def test_filter_sample_ids_exclude(self):
        obs = _filter_contigs(
            contigs=self.contigs,
            on="sample",
            ids=["sample2", "sample3"],
            exclude_ids=True,
        )

        self.assertDictEqual(
            obs.sample_dict(), {"sample1": os.path.join(obs.path, "sample1_contigs.fa")}
        )

    def test_filter_sample_ids_and_metadata_are_combined(self):
        obs = _filter_contigs(
            contigs=self.contigs,
            on="sample",
            ids=["sample2"],
            metadata=self.metadata,
            where="col1='yes'",
        )

        self.assertDictEqual(
            obs.sample_dict(),
            {
                "sample1": os.path.join(obs.path, "sample1_contigs.fa"),
                "sample2": os.path.join(obs.path, "sample2_contigs.fa"),
                "sample3": os.path.join(obs.path, "sample3_contigs.fa"),
            },
        )

    def test_filter_invalid_ids_valid_metadata(self):
        # metadata IDs (sample1, sample3) are present, but `ids` are not
        with self.assertRaisesRegex(
            ValueError, "not present in the contig data: sample4"
        ):
            _filter_contigs(
                contigs=self.contigs,
                on="sample",
                ids=["sample4"],
                metadata=self.metadata,
                where="col1='yes'",
            )

    def test_filter_ids_and_remove_empty(self):
        with self.assertRaisesRegex(ValueError, "No samples remain after filtering"):
            _filter_contigs(
                contigs=self.contigs,
                on="sample",
                ids=["sample3"],
                remove_empty=True,
            )

    def test_filter_metadata_with_extra_ids(self):
        # metadata may list superset of ids
        metadata = q2.Metadata(
            pd.DataFrame(
                data={"col1": ["yes", "yes"]},
                index=pd.Index(["sample1", "sample4"], name="id"),
            )
        )
        obs = _filter_contigs(contigs=self.contigs, on="sample", metadata=metadata)

        self.assertDictEqual(
            obs.sample_dict(), {"sample1": os.path.join(obs.path, "sample1_contigs.fa")}
        )

    def test_filter_contig_id_metadata_on_sample(self):
        metadata = q2.Metadata(
            pd.DataFrame(
                data={"length": [307, 350]},
                index=pd.Index(["k141_0", "k141_2"], name="id"),
            )
        )
        with self.assertRaisesRegex(ValueError, "No samples remain after filtering"):
            _filter_contigs(contigs=self.contigs, on="sample", metadata=metadata)


class TestFilterContigsByContig(FilterContigsTestBase):

    def setUp(self):
        super().setUp()
        self.contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        self.metadata = q2.Metadata(
            pd.DataFrame(
                data={"length": [307, 350, 234]},
                index=pd.Index(["k141_0", "k141_2", "k145_4"], name="id"),
            )
        )
        self.contigs_renamed = ContigSequencesDirFmt(
            self.get_data_path("contigs-renamed"), "r"
        )
        self.metadata_renamed = q2.Metadata(
            pd.DataFrame(
                data={"length": [309, 316, 350]},
                index=pd.Index(
                    [
                        "sample1:dPj72SpWBxMXAVyPtYqJG3",
                        "sample2:2Ta5AWdPGsJiCUZJ4NeP94",
                        "sample2:2Ta5AWdPGsJiCUZJ4NeP92",
                    ],
                    name="id",
                ),
            )
        )
        self.contigs_pooled = ContigSequencesDirFmt(self.get_data_path("pooled"), "r")
        self.contigs_coassembled = ContigSequencesDirFmt(
            self.get_data_path("coassembled"), "r"
        )

    def test_filter_keeps_file_names(self):
        obs = _filter_contigs(contigs=self.contigs, on="contig", ids=["k141_0"])
        self.assertDictEqual(
            obs.sample_dict(),
            {
                "sample1": os.path.join(obs.path, "sample1_contigs.fa"),
                "sample2": os.path.join(obs.path, "sample2_contigs.fa"),
            },
        )

    def test_filter_ids(self):
        obs = _filter_contigs(
            contigs=self.contigs, on="contig", ids=["k141_0", "k145_4"]
        )
        self.assertListEqual(self._read_all_contig_ids(obs), ["k141_0", "k145_4"])

    def test_filter_ids_exclude(self):
        obs = _filter_contigs(
            contigs=self.contigs,
            on="contig",
            ids=["k141_0", "k141_1", "k141_2", "k141_3", "k141_4"],
            exclude_ids=True,
        )
        self.assertListEqual(
            self._read_all_contig_ids(obs),
            [
                "k141_9",
                "k141_6",
                "k141_7",
                "k141_8",
                "k141_5",
                "k145_4",
                "k145_6",
                "k145_7",
                "k145_8",
            ],
        )

    def test_filter_metadata_where(self):
        obs = _filter_contigs(
            contigs=self.contigs,
            on="contig",
            metadata=self.metadata,
            where="length<310",
        )
        self.assertDictEqual(
            obs.sample_dict(),
            {
                "sample1": os.path.join(obs.path, "sample1_contigs.fa"),
                "sample2": os.path.join(obs.path, "sample2_contigs.fa"),
            },
        )
        self.assertListEqual(self._read_all_contig_ids(obs), ["k141_0", "k145_4"])

    def test_filter_ids_and_metadata_are_combined(self):
        obs = _filter_contigs(
            contigs=self.contigs,
            on="contig",
            ids=["k141_9"],
            metadata=self.metadata,
            where="length>300",
        )
        self.assertListEqual(
            self._read_all_contig_ids(obs), ["k141_9", "k141_2", "k141_0"]
        )

    def test_filter_length_and_ids(self):
        obs = _filter_contigs(
            contigs=self.contigs,
            on="contig",
            ids=["k141_0", "k141_2", "k145_4"],
            length_threshold=320,
            remove_empty=True,
        )
        self.assertListEqual(self._read_all_contig_ids(obs), ["k141_2"])
        self.assertDictEqual(
            obs.sample_dict(),
            {
                "sample1": os.path.join(obs.path, "sample1_contigs.fa"),
            },
        )

    def test_filter_length_and_ids_exclude(self):
        obs = _filter_contigs(
            contigs=self.contigs,
            on="contig",
            ids=["k141_2", "k141_8"],
            exclude_ids=True,
            length_threshold=320,
        )
        self.assertListEqual(
            self._read_all_contig_ids(obs),
            ["k141_7", "k141_5", "k141_3", "k141_1", "k145_6"],
        )

    def test_filter_missing_ids(self):
        with self.assertRaisesRegex(
            ValueError, "not present in the contig data: id1, id2"
        ):
            _filter_contigs(
                contigs=self.contigs, on="contig", ids=["k141_0", "id2", "id1"]
            )

    def test_filter_metadata_with_sample_ids(self):
        metadata = q2.Metadata(
            pd.DataFrame(data={"col1": ["yes"]}, index=pd.Index(["sample1"], name="id"))
        )
        with self.assertRaisesRegex(ValueError, "No contigs remain after filtering"):
            _filter_contigs(
                contigs=self.contigs,
                on="contig",
                metadata=metadata,
                remove_empty=True,
            )

    def test_filter_metadata_with_extra_ids(self):
        metadata = q2.Metadata(
            pd.DataFrame(
                data={"length": [307, 1]},
                index=pd.Index(["k141_0", "nope"], name="id"),
            )
        )
        obs = _filter_contigs(contigs=self.contigs, on="contig", metadata=metadata)
        self.assertListEqual(self._read_all_contig_ids(obs), ["k141_0"])

    def test_filter_everything(self):
        with self.assertRaisesRegex(ValueError, "No contigs remain after filtering"):
            _filter_contigs(
                contigs=self.contigs,
                on="contig",
                length_threshold=1000,
                remove_empty=True,
            )

    # renamed SampleData[Contigs]: '<sample_id>:<shortuuid>' in per-sample files
    def test_filter_renamed_whole_and_short_ids(self):
        obs = _filter_contigs(
            contigs=self.contigs_renamed,
            on="contig",
            ids=["sample1:dPj72SpWBxMXAVyPtYqJG1", "2Ta5AWdPGsJiCUZJ4NeP92"],
        )
        self.assertDictEqual(
            {
                "sample1": self._read_contig_ids(obs.sample_dict()["sample1"]),
                "sample2": self._read_contig_ids(obs.sample_dict()["sample2"]),
            },
            {
                "sample1": ["sample1:dPj72SpWBxMXAVyPtYqJG1"],
                "sample2": ["sample2:2Ta5AWdPGsJiCUZJ4NeP92"],
            },
        )

    def test_filter_renamed_metadata_where(self):
        obs = _filter_contigs(
            contigs=self.contigs_renamed,
            on="contig",
            metadata=self.metadata_renamed,
            where="length>310",
            remove_empty=True,
        )
        self.assertDictEqual(
            obs.sample_dict(), {"sample2": os.path.join(obs.path, "sample2_contigs.fa")}
        )
        self.assertListEqual(
            self._read_all_contig_ids(obs),
            ["sample2:2Ta5AWdPGsJiCUZJ4NeP92", "sample2:2Ta5AWdPGsJiCUZJ4NeP94"],
        )

    def test_filter_renamed_wrong_sample_whole_id(self):
        # the shortuuid belongs to a sample1 contig
        with self.assertRaisesRegex(
            ValueError, "not present in the contig data: sample2:dPj72SpWBxMXAVyPtYqJG1"
        ):
            _filter_contigs(
                contigs=self.contigs_renamed,
                on="contig",
                ids=["sample2:dPj72SpWBxMXAVyPtYqJG1"],
            )

    # pooled FeatureData[Contig]: '<sample_id>:<shortuuid>' in a single file
    def test_filter_pooled_whole_and_short_ids(self):
        obs = _filter_contigs(
            contigs=self.contigs_pooled,
            on="contig",
            ids=["sample1:dPj72SpWBxMXAVyPtYqJG1", "2Ta5AWdPGsJiCUZJ4NeP92"],
        )
        self.assertDictEqual(
            obs.sample_dict(), {"pooled": os.path.join(obs.path, "pooled_contigs.fa")}
        )
        self.assertListEqual(
            self._read_all_contig_ids(obs),
            ["sample1:dPj72SpWBxMXAVyPtYqJG1", "sample2:2Ta5AWdPGsJiCUZJ4NeP92"],
        )

    def test_filter_pooled_short_ids_exclude(self):
        obs = _filter_contigs(
            contigs=self.contigs_pooled,
            on="contig",
            ids=[
                "2Ta5AWdPGsJiCUZJ4NeP91",
                "2Ta5AWdPGsJiCUZJ4NeP92",
                "2Ta5AWdPGsJiCUZJ4NeP93",
                "2Ta5AWdPGsJiCUZJ4NeP94",
            ],
            exclude_ids=True,
        )
        self.assertListEqual(self._read_all_contig_ids(obs), POOLED_SAMPLE1_IDS)

    # coassembled FeatureData[Contig]: assembler IDs in a single file
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


class TestFilterContigsBySamplePrefix(FilterContigsTestBase):

    def setUp(self):
        super().setUp()
        self.contigs = ContigSequencesDirFmt(self.get_data_path("pooled"), "r")

    def test_filter_ids(self):
        obs = _filter_contigs(contigs=self.contigs, on="sample_prefix", ids=["sample1"])
        self.assertListEqual(self._read_all_contig_ids(obs), POOLED_SAMPLE1_IDS)

    def test_filter_ids_exclude(self):
        obs = _filter_contigs(
            contigs=self.contigs,
            on="sample_prefix",
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
            contigs=self.contigs,
            on="sample_prefix",
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
            _filter_contigs(contigs=self.contigs, on="sample_prefix", ids=["sample3"])

    def test_sample_prefix_any_separator(self):
        uuid = "mtjebimcR24S9DZ62TY6Fh"
        for separator in [":", ";", "_", "|", ".", "C"]:
            self.assertEqual(
                _match_sample_prefix_ids(f"sample1{separator}{uuid}", {"sample1"}),
                {"sample1"},
            )

    def test_sample_prefix_separator_in_sample_id(self):
        self.assertEqual(
            _match_sample_prefix_ids(
                "my_sample_1_mtjebimcR24S9DZ62TY6Fh", {"my_sample_1"}
            ),
            {"my_sample_1"},
        )

    def test_sample_prefix_not_matched_by_shorter_sample_id(self):
        self.assertEqual(
            _match_sample_prefix_ids("sample10:mtjebimcR24S9DZ62TY6Fh", {"sample1"}),
            set(),
        )

    def test_contig_matched_by_whole_id_or_shortuuid(self):
        contig_id = "sample1:mtjebimcR24S9DZ62TY6Fh"
        self.assertEqual(_match_contig_ids(contig_id, {contig_id}), {contig_id})
        self.assertEqual(
            _match_contig_ids(contig_id, {"mtjebimcR24S9DZ62TY6Fh"}),
            {"mtjebimcR24S9DZ62TY6Fh"},
        )

    def test_contig_not_matched_by_partial_suffix(self):
        self.assertEqual(_match_contig_ids("k141_19", {"9", "k141_9"}), set())

    def test_contig_matches_both_whole_id_and_shortuuid(self):
        contig_id = "sample1:mtjebimcR24S9DZ62TY6Fh"
        selected_ids = {contig_id, "mtjebimcR24S9DZ62TY6Fh", "sample1", "unrelated"}
        self.assertEqual(
            _match_contig_ids(contig_id, selected_ids),
            {contig_id, "mtjebimcR24S9DZ62TY6Fh"},
        )

    def test_filter_metadata_with_contig_ids(self):
        metadata = q2.Metadata(
            pd.DataFrame(
                data={"col1": ["yes"]},
                index=pd.Index(["sample1:dPj72SpWBxMXAVyPtYqJG1"], name="id"),
            )
        )
        with self.assertWarnsRegex(UserWarning, "No contigs remain after filtering"):
            obs = _filter_contigs(
                contigs=self.contigs, on="sample_prefix", metadata=metadata
            )
        self.assertListEqual(self._read_all_contig_ids(obs), [])


class TestFilterContigsPipeline(FilterContigsTestBase):

    def setUp(self):
        super().setUp()
        self.filter_contigs = self.plugin.pipelines["filter_contigs"]
        self.sample_contigs = q2.Artifact.import_data(
            "SampleData[Contigs]", self.get_data_path("contigs-with-empty")
        )
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

    def test_filter_metadata_no_query_no_metadata(self):
        with self.assertRaisesRegex(
            ValueError, "At least one of the following parameters must be provided"
        ):
            self.filter_contigs(contigs=self.sample_contigs)

    def test_empty_metadata_query_include(self):
        with self.assertRaisesRegex(ValueError, "No samples remain after filtering"):
            self.filter_contigs(
                contigs=self.sample_contigs,
                metadata=self.sample_metadata,
                where="selected='no'",
            )

    def test_empty_metadata_query_include_pooled(self):
        with self.assertWarnsRegex(UserWarning, "No contigs remain after filtering"):
            (obs,) = self.filter_contigs(
                contigs=self.pooled_contigs,
                metadata=self.sample_metadata,
                where="selected='no'",
            )
        self.assertDictEqual(
            obs.view(ContigSequencesDirFmt).sample_dict(relative=True),
            self.pooled_contigs.view(ContigSequencesDirFmt).sample_dict(relative=True),
        )
        self.assertListEqual(self._read_all_contig_ids(obs), [])

    def test_empty_metadata_query_exclude(self):
        (obs,) = self.filter_contigs(
            contigs=self.sample_contigs,
            metadata=self.sample_metadata,
            where="selected='no'",
            exclude_ids=True,
        )
        self.assertDictEqual(
            obs.view(ContigSequencesDirFmt).sample_dict(relative=True),
            self.sample_contigs.view(ContigSequencesDirFmt).sample_dict(relative=True),
        )
        self.assertListEqual(
            self._read_all_contig_ids(obs),
            self._read_all_contig_ids(self.sample_contigs),
        )

    def test_empty_metadata_query_combined_with_explicit_ids(self):
        (obs,) = self.filter_contigs(
            contigs=self.sample_contigs,
            ids=["sample1"],
            metadata=self.sample_metadata,
            where="selected='no'",
        )
        self.assertEqual(
            set(obs.view(ContigSequencesDirFmt).sample_dict()), {"sample1"}
        )

    def test_coassembled_contigs_remove_empty(self):
        with self.assertRaisesRegex(ValueError, "No solution for inputs"):
            self.filter_contigs(
                contigs=self.coassembled_contigs,
                on="contig",
                length_threshold=300,
                remove_empty=True,
            )

    def test_pooled_contigs_remove_empty(self):
        # pooled contigs are a single file - only remove_empty=False is allowed
        with self.assertRaisesRegex(ValueError, "No solution for inputs"):
            self.filter_contigs(
                contigs=self.pooled_contigs, length_threshold=300, remove_empty=True
            )

    def test_coassembled_contigs_requires_explicit_on(self):
        with self.assertRaisesRegex(ValueError, "No solution for inputs"):
            self.filter_contigs(contigs=self.coassembled_contigs, length_threshold=300)

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

    def test_where_without_metadata(self):
        with self.assertRaisesRegex(ValueError, "Metadata must be provided"):
            self.filter_contigs(
                contigs=self.sample_contigs,
                length_threshold=100,
                where="x='y'",
            )

    def test_exclude_without_ids(self):
        with self.assertRaisesRegex(ValueError, "Either 'ids' or metadata"):
            self.filter_contigs(
                contigs=self.sample_contigs,
                length_threshold=100,
                exclude_ids=True,
            )
