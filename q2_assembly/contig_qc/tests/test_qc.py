import os
from pathlib import Path

import numpy as np
import pandas as pd
import qiime2 as q2
from q2_types.per_sample_sequences import ContigSequencesDirFmt
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import assembly

from q2_assembly.contig_qc.qc import (
    _process_single_fasta,
    _get_gc_content,
    _get_contig_lengths,
    _calculate_cumulative_length,
    _calculate_nx_metrics,
    _calculate_all_metrics,
    compute_sample_metrics,
    dump_all_to_arrow,
    _reset_indices,
    process_metadata,
    _evaluate_contigs,
)


class TestCoreCalculations(TestPluginBase):
    package = "q2_assembly.contig_qc.tests"

    def setUp(self):
        super().setUp()
        self.raw_data_sampleD = {
            "lengths": [4, 2, 5],
            "total_length": 11,
            "num_contigs": 3,
            "gc": [50.0, 100.0, 80.0],
            "sorted_lengths": [5, 4, 2],
        }
        self.raw_data_empty_sample = {
            "lengths": [],
            "total_length": 0,
            "num_contigs": 0,
            "gc": [],
            "sorted_lengths": [],
        }

    def test_process_single_fasta_ok(self):
        fp = Path(self.get_data_path("sampleA_ok.fa"))
        sample_id, data = _process_single_fasta(fp)
        self.assertEqual(sample_id, "sampleA_ok")
        expected_data = {
            "lengths": [4, 2, 5],
            "total_length": 11,
            "num_contigs": 3,
            "gc": [50.0, 100.0, 80.0],
            "sorted_lengths": [5, 4, 2],
        }
        self.assertDictEqual(data, expected_data)

    def test_process_single_fasta_empty_file(self):
        sample_id, data = _process_single_fasta(
            Path(self.get_data_path("sampleB_empty.fa"))
        )
        self.assertEqual(sample_id, "sampleB_empty")
        expected_data = {
            "lengths": [],
            "total_length": 0,
            "num_contigs": 0,
            "gc": [],
            "sorted_lengths": [],
        }
        self.assertDictEqual(data, expected_data)

    def test_get_gc_content_ok(self):
        rows = _get_gc_content("sampleD", self.raw_data_sampleD)
        expected_rows = [
            {"sample": "sampleD", "gc": 50.0},
            {"sample": "sampleD", "gc": 100.0},
            {"sample": "sampleD", "gc": 80.0},
        ]
        self.assertEqual(len(rows), len(expected_rows))
        for obs, exp in zip(rows, expected_rows):
            self.assertDictEqual(obs, exp)

    def test_get_gc_content_empty(self):
        rows = _get_gc_content("sample_empty", self.raw_data_empty_sample)
        self.assertEqual(rows, [])

    def test_get_contig_lengths_ok(self):
        rows = _get_contig_lengths("sampleD", self.raw_data_sampleD)
        expected_rows = [
            {"sample": "sampleD", "contig_length": 4},
            {"sample": "sampleD", "contig_length": 2},
            {"sample": "sampleD", "contig_length": 5},
        ]
        self.assertEqual(len(rows), len(expected_rows))
        for obs, exp in zip(rows, expected_rows):
            self.assertDictEqual(obs, exp)

    def test_get_contig_lengths_empty(self):
        rows = _get_contig_lengths("sample_empty", self.raw_data_empty_sample)
        self.assertEqual(rows, [])

    def test_calculate_cumulative_length_no_pruning(self):
        sorted_lengths = [100, 80, 50, 20]
        df = _calculate_cumulative_length("s1", sorted_lengths)
        self.assertEqual(len(df), 4)
        expected_df = pd.DataFrame(
            {
                "sample": ["s1"] * 4,
                "rank": [1, 2, 3, 4],
                "cumulative_length": [100, 180, 230, 250],
            }
        )
        pd.testing.assert_frame_equal(df, expected_df)

    def test_calculate_cumulative_length_with_pruning(self):
        sorted_lengths = list(range(600, 0, -1))
        df = _calculate_cumulative_length("s2", sorted_lengths, max_points=10)
        self.assertEqual(len(df), 10)
        self.assertEqual(df["rank"].iloc[0], 1)
        self.assertEqual(df["rank"].iloc[-1], 600)
        self.assertEqual(df["cumulative_length"].iloc[0], sorted_lengths[0])
        self.assertEqual(df["cumulative_length"].iloc[-1], sum(sorted_lengths))

    def test_calculate_cumulative_length_empty(self):
        df = _calculate_cumulative_length("s3", [])
        self.assertIsNone(df)

    def test_calculate_cumulative_length_exactly_max_points(self):
        sorted_lengths = list(range(10, 0, -1))
        df = _calculate_cumulative_length("s4", sorted_lengths, max_points=10)
        self.assertEqual(len(df), 10)
        self.assertEqual(df["rank"].iloc[-1], 10)
        self.assertEqual(df["cumulative_length"].iloc[-1], sum(sorted_lengths))

    def test_calculate_nx_metrics_ok(self):
        sorted_lengths = [100, 80, 60, 40, 20]
        total_length = 300
        nx_rows = _calculate_nx_metrics("s1", sorted_lengths, total_length)
        self.assertEqual(len(nx_rows), 100)
        self.assertIn({"sample": "s1", "percent": 50, "nx": 80, "lx": 2}, nx_rows)
        self.assertIn({"sample": "s1", "percent": 90, "nx": 40, "lx": 4}, nx_rows)
        self.assertIn({"sample": "s1", "percent": 100, "nx": 20, "lx": 5}, nx_rows)

    def test_calculate_nx_metrics_empty_lengths(self):
        nx_rows = _calculate_nx_metrics("s2", [], 0)
        self.assertEqual(nx_rows, [])

    def test_calculate_nx_metrics_total_zero(self):
        nx_rows = _calculate_nx_metrics("s3", [10, 5], 0)
        self.assertEqual(nx_rows, [])

    def test_calculate_all_metrics_ok(self):
        sample_id = "sampleD"
        all_metrics = _calculate_all_metrics(sample_id, self.raw_data_sampleD)

        expected_gc_rows = [
            {"sample": "sampleD", "gc": 50.0},
            {"sample": "sampleD", "gc": 100.0},
            {"sample": "sampleD", "gc": 80.0},
        ]
        expected_len_rows = [
            {"sample": "sampleD", "contig_length": 4},
            {"sample": "sampleD", "contig_length": 2},
            {"sample": "sampleD", "contig_length": 5},
        ]
        expected_cumul_df = pd.DataFrame(
            {
                "sample": ["sampleD"] * 3,
                "rank": [1, 2, 3],
                "cumulative_length": [5, 9, 11],
            }
        )

        self.assertEqual(len(all_metrics["seq_gc_rows"]), len(expected_gc_rows))
        for obs, exp in zip(all_metrics["seq_gc_rows"], expected_gc_rows):
            self.assertDictEqual(obs, exp)

        self.assertEqual(len(all_metrics["seq_len_rows"]), len(expected_len_rows))
        for obs, exp in zip(all_metrics["seq_len_rows"], expected_len_rows):
            self.assertDictEqual(obs, exp)

        pd.testing.assert_frame_equal(
            all_metrics["cumulative_df_part"], expected_cumul_df
        )

        self.assertIs(len(all_metrics["nx_rows"]), 100)

    def test_calculate_all_metrics_empty(self):
        sample_id = "sample_empty"
        all_metrics = _calculate_all_metrics(sample_id, self.raw_data_empty_sample)
        expected_metrics = {
            "seq_gc_rows": [],
            "seq_len_rows": [],
            "cumulative_df_part": None,
            "nx_rows": [],
        }
        self.assertDictEqual(all_metrics, expected_metrics)


class TestSummariesAndHelpers(TestPluginBase):
    package = "q2_assembly.contig_qc.tests"

    def setUp(self):
        super().setUp()
        self.sample_contigs_df = pd.read_csv(
            self.get_data_path("sample_contigs_df.csv")
        )

    def test_compute_sample_metrics(self):
        observed_metrics = compute_sample_metrics(self.sample_contigs_df)
        expected_metrics = pd.DataFrame.from_records(
            [
                {
                    "sample": "s1",
                    "count": 3,
                    "mean": np.round((100 + 50 + 25) / 3),
                    "longest": 100,
                    "n50": 100,
                    "n90": 25,
                    "l50": 1,
                    "l90": 3,
                    "total_length": 175,
                },
                {
                    "sample": "s2",
                    "count": 2,
                    "mean": 150.0,
                    "longest": 200,
                    "n50": 200,
                    "n90": 100,
                    "l50": 1,
                    "l90": 2,
                    "total_length": 300,
                },
            ]
        ).set_index("sample", drop=False)
        expected_metrics.index.name = "id"
        pd.testing.assert_frame_equal(observed_metrics, expected_metrics)

    def test_compute_sample_metrics_empty_df(self):
        empty_df = pd.DataFrame(columns=["sample", "contig_length", "meta_A"])
        observed_metrics = compute_sample_metrics(empty_df)
        pd.testing.assert_frame_equal(
            observed_metrics,
            pd.DataFrame(
                columns=[
                    "sample",
                    "count",
                    "mean",
                    "longest",
                    "n50",
                    "n90",
                    "l50",
                    "l90",
                    "total_length",
                ],
                index=pd.Index([], name="id"),
            ),
        )

    def test_compute_sample_metrics_single_contig(self):
        single_contig_df = pd.DataFrame(
            {"sample": ["s3"], "contig_length": [500], "meta_C": ["z"]}
        )
        observed_metrics = compute_sample_metrics(single_contig_df)
        expected_metrics = pd.DataFrame.from_records(
            [
                {
                    "sample": "s3",
                    "count": 1,
                    "mean": 500.0,
                    "longest": 500,
                    "n50": 500,
                    "n90": 500,
                    "l50": 1,
                    "l90": 1,
                    "total_length": 500,
                }
            ]
        ).set_index("sample", drop=False)
        expected_metrics.index.name = "id"
        pd.testing.assert_frame_equal(observed_metrics, expected_metrics)

    def test_dump_all_to_arrow(self):
        data = {
            "seq_len_df": pd.DataFrame(),
            "seq_gc_df": pd.DataFrame(),
            "cumulative_df": pd.DataFrame(),
            "nx_df": pd.DataFrame(),
        }
        dump_all_to_arrow(data, self.temp_dir.name)
        for fn in [
            "contig_length_data.arrow",
            "gc_content_data.arrow",
            "cumulative_length_data.arrow",
            "nx_curve_data.arrow",
        ]:
            self.assertTrue(Path(os.path.join(self.temp_dir.name, "data", fn)).exists())

    def test_reset_indices(self):
        df1 = q2.Artifact.import_data(
            "ImmutableMetadata",
            q2.Metadata(
                pd.DataFrame(
                    {"col1": ["a", "b"]}, index=pd.Index(["s1", "s2"], name="sample-id")
                )
            ),
        )
        df2 = q2.Artifact.import_data(
            "ImmutableMetadata",
            q2.Metadata(
                pd.DataFrame(
                    {"col1": ["c", "d"]}, index=pd.Index(["s1", "s2"], name="sample-id")
                )
            ),
        )
        observed_df = _reset_indices([df1, df2]).view(q2.Metadata).to_dataframe()
        expected_df = pd.DataFrame(
            {"col1": ["a", "b", "c", "d"]},
            index=pd.Index(["0", "1", "2", "3"], name="id"),
        )
        pd.testing.assert_frame_equal(observed_df, expected_df)

    def test_process_metadata(self):
        meta = q2.Metadata(
            pd.read_csv(self.get_data_path("metadata.tsv"), sep="\t", index_col=0)
        )
        samples_by_metadata = {"all_samples": ["s1", "s2", "s3"]}

        obs_values, obs_sample_dict = process_metadata(meta, samples_by_metadata)

        self.assertDictEqual(obs_values, {"cat1": ["a", "b"], "cat2": ["abc", "def"]})
        self.assertDictEqual(
            obs_sample_dict,
            {
                "all_samples": ["s1", "s2", "s3"],
                "cat1": {"a": ["s1", "s3"], "b": ["s2"]},
                "cat2": {"abc": ["s1", "s2"], "def": ["s3"]},
            },
        )


class TestIntegration(TestPluginBase):
    package = "q2_assembly.contig_qc.tests"

    def setUp(self):
        super().setUp()
        self.contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        self.contigs_artifact = q2.Artifact.import_data(
            "SampleData[Contigs]",
            self.contigs,
        )
        self.metadata = q2.Metadata(
            pd.read_csv(self.get_data_path("metadata_small.tsv"), sep="\t", index_col=0)
        )

    def test_evaluate_contigs_helper(self):
        obs_results = _evaluate_contigs(self.contigs, n_cpus=1)
        exp_results = (
            pd.read_csv(self.get_data_path("sample_metrics.csv"), index_col=0),
            pd.read_csv(self.get_data_path("nx.csv"), index_col=0),
            pd.read_csv(self.get_data_path("gc.csv"), index_col=0),
            pd.read_csv(self.get_data_path("length.csv"), index_col=0),
            pd.read_csv(self.get_data_path("cumulative.csv"), index_col=0),
        )

        for obs, exp in zip(obs_results, exp_results):
            exp.index = exp.index.astype(str)

        obs_metrics, obs_nx, obs_gc, obs_len, obs_cumul = obs_results
        exp_metrics, exp_nx, exp_gc, exp_len, exp_cumul = exp_results

        # test sample metrics
        pd.testing.assert_frame_equal(
            obs_metrics.to_dataframe(), exp_metrics, check_dtype=False
        )

        # test Nx metrics
        obs_nx = obs_nx.to_dataframe()
        obs_nx.sort_values(by=["sample", "percent"], inplace=True).reset_index(
            inplace=True
        )
        pd.testing.assert_frame_equal(obs_nx, exp_nx, check_dtype=False)

        # test GC metrics
        pd.testing.assert_frame_equal(obs_gc.to_dataframe(), exp_gc, check_dtype=False)

        # test contig lengths
        pd.testing.assert_frame_equal(
            obs_len.to_dataframe(), exp_len, check_dtype=False
        )

        # test cumulative lengths
        pd.testing.assert_frame_equal(
            obs_cumul.to_dataframe(), exp_cumul, check_dtype=False
        )

    def test_evaluate_contigs_pipeline_single_partition(self):
        obs_results, obs_viz = assembly.pipelines.evaluate_contigs(
            contigs=self.contigs_artifact,
            metadata=self.metadata,
            num_partitions=1,
        )

        # assert the results table
        self.assertIsInstance(obs_results, q2.Artifact)
        exp_results = pd.read_csv(self.get_data_path("results_table.csv"), index_col=0)
        pd.testing.assert_frame_equal(
            obs_results.view(q2.Metadata).to_dataframe(), exp_results
        )

        # assert the visualization
        self.assertIsInstance(obs_viz, q2.Visualization)
        obs_viz.export_data(self.temp_dir.name)
        for page in ("index.html", "table.html", "grouped.html"):
            self.assertTrue(Path(os.path.join(self.temp_dir.name, page)).exists())
        with open(
            os.path.join(self.temp_dir.name, "index.html"), "r", encoding="utf-8"
        ) as f:
            html_content = f.read()
            self.assertIn('id="metadata-card"', html_content)
        for arrow in (
            "contig_length_data.arrow",
            "gc_content_data.arrow",
            "cumulative_length_data.arrow",
            "nx_curve_data.arrow",
        ):
            self.assertTrue(
                Path(os.path.join(self.temp_dir.name, "data", arrow)).exists()
            )
        pd.testing.assert_frame_equal(
            pd.read_csv(
                os.path.join(self.temp_dir.name, "data", "sample_metrics.tsv"),
                sep="\t",
                index_col=0,
            ),
            exp_results.merge(
                self.metadata.filter_columns(column_type="categorical").to_dataframe(),
                left_on="sample",
                right_index=True,
            ),
        )

    def test_evaluate_contigs_pipeline_two_partitions(self):
        obs_results, obs_viz = assembly.pipelines.evaluate_contigs(
            contigs=self.contigs_artifact,
            metadata=self.metadata,
            num_partitions=2,
        )

        self.assertIsInstance(obs_results, q2.Artifact)
        self.assertIsInstance(obs_viz, q2.Visualization)

        exp_results = pd.read_csv(self.get_data_path("results_table.csv"), index_col=0)
        pd.testing.assert_frame_equal(
            obs_results.view(q2.Metadata).to_dataframe(), exp_results
        )

    def test_evaluate_contigs_pipeline_no_metadata(self):
        obs_results, obs_viz = assembly.pipelines.evaluate_contigs(
            contigs=self.contigs_artifact,
            num_partitions=1,
        )

        # assert the results table
        self.assertIsInstance(obs_results, q2.Artifact)
        exp_results = pd.read_csv(self.get_data_path("results_table.csv"), index_col=0)
        pd.testing.assert_frame_equal(
            obs_results.view(q2.Metadata).to_dataframe(), exp_results
        )

        # assert the visualization
        self.assertIsInstance(obs_viz, q2.Visualization)
        obs_viz.export_data(self.temp_dir.name)
        for page in ("index.html", "table.html"):
            self.assertTrue(Path(os.path.join(self.temp_dir.name, page)).exists())
        self.assertFalse(
            Path(os.path.join(self.temp_dir.name, "grouped.html")).exists()
        )
        with open(
            os.path.join(self.temp_dir.name, "index.html"), "r", encoding="utf-8"
        ) as f:
            html_content = f.read()
            self.assertNotIn('id="metadata-card"', html_content)
