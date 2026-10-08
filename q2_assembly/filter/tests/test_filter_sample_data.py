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

from q2_assembly.filter.filter import _filter_contigs
from q2_assembly.filter.tests._utils import FilterContigsTestBase


class TestFilterContigsBySample(FilterContigsTestBase):

    def setUp(self):
        super().setUp()
        self.sample_data_contigs = ContigSequencesDirFmt(
            self.get_data_path("contigs-with-empty"), "r"
        )
        self.metadata_df = pd.DataFrame(
            data={"col1": ["yes", "no", "yes"]},
            index=pd.Index(["sample1", "sample2", "sample3"], name="id"),
        )
        self.metadata = q2.Metadata(self.metadata_df)

    def test_filter_sample_metadata(self):
        obs = _filter_contigs(
            contigs=self.sample_data_contigs,
            on="sample",
            filter_samples=True,
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
            contigs=self.sample_data_contigs,
            on="sample",
            filter_samples=True,
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
        obs = _filter_contigs(
            contigs=self.sample_data_contigs,
            on="sample",
            filter_samples=True,
            metadata=self.metadata,
        )

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
            contigs=self.sample_data_contigs,
            on="sample",
            filter_samples=True,
            metadata=metadata_subset,
        )

        self.assertDictEqual(
            obs.sample_dict(),
            {
                "sample1": os.path.join(obs.path, "sample1_contigs.fa"),
                "sample2": os.path.join(obs.path, "sample2_contigs.fa"),
            },
        )

    def test_filter_by_length(self):
        obs = _filter_contigs(
            contigs=self.sample_data_contigs,
            on="sample",
            filter_samples=True,
            length_threshold=320,
        )

        obs_counts = {
            _id: len(self._read_contig_ids(fp)) for _id, fp in obs.sample_dict().items()
        }
        self.assertDictEqual(obs_counts, {"sample1": 6, "sample2": 1, "sample3": 0})

    def test_filter_remove_empty(self):
        obs = _filter_contigs(
            contigs=self.sample_data_contigs,
            on="sample",
            filter_samples=True,
            remove_empty=True,
        )

        self.assertDictEqual(
            obs.sample_dict(),
            {
                "sample1": os.path.join(obs.path, "sample1_contigs.fa"),
                "sample2": os.path.join(obs.path, "sample2_contigs.fa"),
            },
        )

    def test_filter_by_length_and_remove_empty(self):
        obs = _filter_contigs(
            contigs=self.sample_data_contigs,
            on="sample",
            filter_samples=True,
            length_threshold=400,
            remove_empty=True,
        )

        self.assertEqual(len(obs.sample_dict()), 1)

        with open(os.path.join(obs.path, "sample1_contigs.fa")) as f:
            self.assertEqual(len(list(skbio.io.read(f, format="fasta"))), 2)

    def test_filter_everything(self):
        with self.assertRaisesRegex(ValueError, "No samples remain after filtering"):
            _filter_contigs(
                contigs=self.sample_data_contigs,
                on="sample",
                filter_samples=True,
                length_threshold=1000,
                remove_empty=True,
            )

    def test_filter_sample_ids(self):
        obs = _filter_contigs(
            contigs=self.sample_data_contigs,
            on="sample",
            filter_samples=True,
            ids=["sample2"],
        )

        self.assertDictEqual(
            obs.sample_dict(), {"sample2": os.path.join(obs.path, "sample2_contigs.fa")}
        )

    def test_filter_sample_ids_exclude(self):
        obs = _filter_contigs(
            contigs=self.sample_data_contigs,
            on="sample",
            filter_samples=True,
            ids=["sample2", "sample3"],
            exclude_ids=True,
        )

        self.assertDictEqual(
            obs.sample_dict(), {"sample1": os.path.join(obs.path, "sample1_contigs.fa")}
        )

    def test_filter_sample_ids_and_metadata_are_combined(self):
        obs = _filter_contigs(
            contigs=self.sample_data_contigs,
            on="sample",
            filter_samples=True,
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
                contigs=self.sample_data_contigs,
                on="sample",
                filter_samples=True,
                ids=["sample4"],
                metadata=self.metadata,
                where="col1='yes'",
            )

    def test_filter_ids_and_remove_empty(self):
        with self.assertRaisesRegex(ValueError, "No samples remain after filtering"):
            _filter_contigs(
                contigs=self.sample_data_contigs,
                on="sample",
                filter_samples=True,
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
        obs = _filter_contigs(
            contigs=self.sample_data_contigs,
            on="sample",
            filter_samples=True,
            metadata=metadata,
        )

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
            _filter_contigs(
                contigs=self.sample_data_contigs,
                on="sample",
                filter_samples=True,
                metadata=metadata,
            )

    def test_empty_metadata_query_include(self):
        with self.assertRaisesRegex(ValueError, "No samples remain after filtering"):
            _filter_contigs(
                contigs=self.sample_data_contigs,
                on="sample",
                filter_samples=True,
                metadata=self.metadata,
                where="col1='missing'",
            )

    def test_empty_metadata_query_exclude(self):
        obs = _filter_contigs(
            contigs=self.sample_data_contigs,
            on="sample",
            filter_samples=True,
            metadata=self.metadata,
            where="col1='missing'",
            exclude_ids=True,
        )
        self.assertDictEqual(
            obs.sample_dict(relative=True),
            self.sample_data_contigs.sample_dict(relative=True),
        )
        self.assertListEqual(
            self._read_all_contig_ids(obs),
            self._read_all_contig_ids(self.sample_data_contigs),
        )

    def test_empty_metadata_query_combined_with_explicit_ids(self):
        obs = _filter_contigs(
            contigs=self.sample_data_contigs,
            on="sample",
            filter_samples=True,
            ids=["sample1"],
            metadata=self.metadata,
            where="col1='missing'",
        )
        self.assertDictEqual(
            obs.sample_dict(), {"sample1": os.path.join(obs.path, "sample1_contigs.fa")}
        )
        self.assertListEqual(
            self._read_all_contig_ids(obs),
            self._read_contig_ids(self.sample_data_contigs.sample_dict()["sample1"]),
        )


class TestFilterSampleDataContigsByContig(FilterContigsTestBase):

    def setUp(self):
        super().setUp()
        self.sample_data_contigs = ContigSequencesDirFmt(
            self.get_data_path("contigs"), "r"
        )
        self.metadata = q2.Metadata(
            pd.DataFrame(
                data={
                    "col1": ["yes", "no", "yes"],
                    "col2": ["yes", "yes", "no"],
                },
                index=pd.Index(["k141_0", "k141_2", "k145_4"], name="id"),
            )
        )
        self.contigs_renamed = ContigSequencesDirFmt(
            self.get_data_path("contigs-renamed"), "r"
        )

    def test_filter_keeps_file_names(self):
        obs = _filter_contigs(
            contigs=self.sample_data_contigs, on="contig", ids=["k141_0"]
        )
        self.assertDictEqual(
            obs.sample_dict(),
            {
                "sample1": os.path.join(obs.path, "sample1_contigs.fa"),
                "sample2": os.path.join(obs.path, "sample2_contigs.fa"),
            },
        )

    def test_filter_ids(self):
        obs = _filter_contigs(
            contigs=self.sample_data_contigs, on="contig", ids=["k141_0", "k145_4"]
        )
        self.assertListEqual(self._read_all_contig_ids(obs), ["k141_0", "k145_4"])

    def test_filter_ids_exclude(self):
        obs = _filter_contigs(
            contigs=self.sample_data_contigs,
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
            contigs=self.sample_data_contigs,
            on="contig",
            metadata=self.metadata,
            where="col1='yes'",
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
            contigs=self.sample_data_contigs,
            on="contig",
            ids=["k141_9"],
            metadata=self.metadata,
            where="col2='yes'",
        )
        self.assertListEqual(
            self._read_all_contig_ids(obs), ["k141_9", "k141_2", "k141_0"]
        )

    def test_filter_length_and_ids(self):
        obs = _filter_contigs(
            contigs=self.sample_data_contigs,
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
            contigs=self.sample_data_contigs,
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
                contigs=self.sample_data_contigs,
                on="contig",
                ids=["k141_0", "id2", "id1"],
            )

    def test_filter_metadata_with_sample_ids(self):
        metadata = q2.Metadata(
            pd.DataFrame(data={"col1": ["yes"]}, index=pd.Index(["sample1"], name="id"))
        )
        with self.assertRaisesRegex(ValueError, "No contigs remain after filtering"):
            _filter_contigs(
                contigs=self.sample_data_contigs,
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
        obs = _filter_contigs(
            contigs=self.sample_data_contigs, on="contig", metadata=metadata
        )
        self.assertListEqual(self._read_all_contig_ids(obs), ["k141_0"])

    def test_filter_everything(self):
        with self.assertRaisesRegex(ValueError, "No contigs remain after filtering"):
            _filter_contigs(
                contigs=self.sample_data_contigs,
                on="contig",
                length_threshold=1000,
                remove_empty=True,
            )

    def test_filter_renamed_whole_and_short_ids(self):
        obs = _filter_contigs(
            contigs=self.contigs_renamed,
            on="contig",
            ids=["sample1:dPj72SpWBxMXAVyPtYqJGA", "2Ta5AWdPGsJiCUZJ4NeP92"],
        )
        self.assertDictEqual(
            {
                "sample1": self._read_contig_ids(obs.sample_dict()["sample1"]),
                "sample2": self._read_contig_ids(obs.sample_dict()["sample2"]),
            },
            {
                "sample1": ["sample1:dPj72SpWBxMXAVyPtYqJGA"],
                "sample2": ["sample2:2Ta5AWdPGsJiCUZJ4NeP92"],
            },
        )

    def test_filter_renamed_wrong_sample_whole_id(self):
        with self.assertRaisesRegex(
            ValueError, "not present in the contig data: sample2:dPj72SpWBxMXAVyPtYqJGA"
        ):
            _filter_contigs(
                contigs=self.contigs_renamed,
                on="contig",
                ids=["sample2:dPj72SpWBxMXAVyPtYqJGA"],
            )


class TestFilterSampleDataContigsPipeline(FilterContigsTestBase):

    def setUp(self):
        super().setUp()
        self.filter_contigs = self.plugin.pipelines["filter_contigs"]
        self.sample_contigs = q2.Artifact.import_data(
            "SampleData[Contigs]", self.get_data_path("contigs-with-empty")
        )
        self.sample_metadata = q2.Metadata(
            pd.DataFrame(
                {"selected": ["yes"]},
                index=pd.Index(["sample1"], name="id"),
            )
        )

    def test_no_filter_params(self):
        with self.assertRaisesRegex(
            ValueError, "At least one of the following parameters must be provided"
        ):
            self.filter_contigs(contigs=self.sample_contigs)

    def test_sample_contigs_on_sample(self):
        (obs,) = self.filter_contigs(
            contigs=self.sample_contigs,
            metadata=self.sample_metadata,
            where="selected='yes'",
        )
        obs_contigs = obs.view(ContigSequencesDirFmt)
        self.assertEqual(obs.type, self.sample_contigs.type)
        self.assertDictEqual(
            obs_contigs.sample_dict(),
            {"sample1": os.path.join(obs_contigs.path, "sample1_contigs.fa")},
        )
        self.assertListEqual(
            self._read_all_contig_ids(obs_contigs),
            self._read_contig_ids(
                self.sample_contigs.view(ContigSequencesDirFmt).sample_dict()["sample1"]
            ),
        )

    def test_sample_contigs_on_contig(self):
        (obs,) = self.filter_contigs(
            contigs=self.sample_contigs,
            on="contig",
            ids=["k141_2", "k145_4"],
        )
        obs_contigs = obs.view(ContigSequencesDirFmt)
        self.assertEqual(obs.type, self.sample_contigs.type)
        self.assertDictEqual(
            obs_contigs.sample_dict(),
            {
                "sample1": os.path.join(obs_contigs.path, "sample1_contigs.fa"),
                "sample2": os.path.join(obs_contigs.path, "sample2_contigs.fa"),
                "sample3": os.path.join(obs_contigs.path, "sample3_contigs.fa"),
            },
        )
        self.assertListEqual(
            self._read_all_contig_ids(obs_contigs), ["k141_2", "k145_4"]
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
