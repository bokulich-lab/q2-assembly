import json
import os
import tempfile
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd
from parameterized import parameterized
from q2_types.per_sample_sequences import ContigSequencesDirFmt
from qiime2 import Metadata
from qiime2.plugin.testing import TestPluginBase

from q2_assembly.contig_qc.qc import (
    _process_single_fasta,
    _get_gc_content,
    _get_contig_lengths,
    _calculate_cumulative_length,
    _calculate_nx_metrics,
    _calculate_all_metrics,
    generate_plotting_data,
    compute_sample_metrics,
    render_spec,
    estimate_column_count,
    evaluate_contigs,
    TEMPLATES,
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

    def test_compute_sample_metrics_no_categories(self):
        observed_metrics = compute_sample_metrics(self.sample_contigs_df, [])
        expected_metrics = [
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
        for obs, exp in zip(
            sorted(observed_metrics, key=lambda x: x["sample"]), expected_metrics
        ):
            self.assertDictEqual(obs, exp)

    def test_compute_sample_metrics_with_categories(self):
        observed_metrics = compute_sample_metrics(
            self.sample_contigs_df, ["meta_A", "meta_B"]
        )
        expected_metrics = [
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
                "meta_A": "x",
                "meta_B": 1,
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
                "meta_A": "y",
                "meta_B": 2,
            },
        ]
        for obs, exp in zip(
            sorted(observed_metrics, key=lambda x: x["sample"]), expected_metrics
        ):
            self.assertDictEqual(obs, exp)

    def test_compute_sample_metrics_empty_df(self):
        empty_df = pd.DataFrame(columns=["sample", "contig_length", "meta_A"])
        metrics = compute_sample_metrics(empty_df, ["meta_A"])
        self.assertEqual(metrics, [])

    def test_compute_sample_metrics_single_contig(self):
        single_contig_df = pd.DataFrame(
            {"sample": ["s3"], "contig_length": [500], "meta_C": ["z"]}
        )
        metrics = compute_sample_metrics(single_contig_df, ["meta_C"])
        self.assertEqual(len(metrics), 1)
        s3_metrics_exp = {
            "sample": "s3",
            "count": 1,
            "mean": 500.0,
            "longest": 500,
            "n50": 500,
            "n90": 500,
            "l50": 1,
            "l90": 1,
            "total_length": 500,
            "meta_C": "z",
        }
        self.assertDictEqual(metrics[0], s3_metrics_exp)

    @parameterized.expand(
        [
            (set(), 4),
            ({"s1", "s2"}, 4),
            ({"sample_long"}, 3),
            ({"sample_very_long_indeed"}, 2),
            ({"short", "another_one_quite_long"}, 2),
        ],
    )
    def test_estimate_column_count(self, sample_ids, expected_cols):
        self.assertEqual(
            estimate_column_count(sample_ids if sample_ids else {""}), expected_cols
        )

    @mock.patch("builtins.open", new_callable=mock.mock_open)
    @mock.patch("q2_assembly.contig_qc.qc.jinja2.Template")
    def test_render_spec(self, mock_template_cls, mock_open_file):
        mock_template_instance = mock.Mock()
        mock_template_instance.render.return_value = "rendered_spec_content"
        mock_template_cls.return_value = mock_template_instance
        expected_path = os.path.join(
            TEMPLATES, "contig_qc", "vega", "my_template.json.j2"
        )
        spec_content = render_spec("my_template.json.j2", key="value")
        mock_open_file.assert_called_once_with(expected_path)
        mock_template_cls.assert_called_once_with(mock_open_file().read())
        mock_template_instance.render.assert_called_once_with(key="value")
        self.assertEqual(spec_content, "rendered_spec_content")


class TestIntegration(TestPluginBase):
    package = "q2_assembly.contig_qc.tests"

    def setUp(self):
        super().setUp()
        self.raw_data_sampleD = {
            "lengths": [10, 20],
            "total_length": 30,
            "num_contigs": 2,
            "gc": [50.0, 50.0],
            "sorted_lengths": [20, 10],
        }

    @mock.patch("q2_assembly.contig_qc.qc.Pool")
    def test_generate_plotting_data_no_metadata(self, mock_pool):
        mock_process_output = [
            ("s1", self.raw_data_sampleD),
            ("s2", self.raw_data_sampleD),
        ]
        mock_calc_output_s1 = _calculate_all_metrics("s1", self.raw_data_sampleD)
        mock_calc_output_s2 = _calculate_all_metrics("s2", self.raw_data_sampleD)
        mock_map = mock.Mock(return_value=mock_process_output)
        mock_starmap = mock.Mock(
            return_value=[mock_calc_output_s1, mock_calc_output_s2]
        )
        mock_pool.return_value.__enter__.return_value.map = mock_map
        mock_pool.return_value.__enter__.return_value.starmap = mock_starmap

        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        results = generate_plotting_data(contigs, metadata_df=None, n_cpus=1)
        self.assertIn("seq_gc_df", results)
        self.assertEqual(
            len(results["seq_gc_df"]), len(mock_calc_output_s1["seq_gc_rows"]) * 2
        )
        self.assertIn("seq_len_df", results)
        self.assertEqual(
            len(results["seq_len_df"]), len(mock_calc_output_s1["seq_len_rows"]) * 2
        )
        self.assertIn("cumulative_df", results)
        self.assertEqual(
            len(results["cumulative_df"]),
            len(mock_calc_output_s1["cumulative_df_part"]) * 2,
        )
        self.assertIn("nx_df", results)
        self.assertEqual(len(results["nx_df"]), len(mock_calc_output_s1["nx_rows"]) * 2)
        self.assertEqual(results["metadata_columns"], [])
        mock_map.assert_called_once()
        mock_starmap.assert_called_once()

    @mock.patch("q2_assembly.contig_qc.qc.Pool")
    def test_generate_plotting_data_with_metadata(self, mock_pool):
        mock_process_output = [("s1", self.raw_data_sampleD)]
        mock_calc_output_s1 = _calculate_all_metrics("s1", self.raw_data_sampleD)
        mock_map = mock.Mock(return_value=mock_process_output)
        mock_starmap = mock.Mock(return_value=[mock_calc_output_s1])
        mock_pool.return_value.__enter__.return_value.map = mock_map
        mock_pool.return_value.__enter__.return_value.starmap = mock_starmap

        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        metadata_df = pd.DataFrame(
            {"meta1": ["val1"], "meta2": ["valA"]},
            index=pd.Index(["s1"], name="sample"),
        )
        results = generate_plotting_data(contigs, metadata_df=metadata_df, n_cpus=1)
        self.assertIn("meta1", results["seq_len_df"].columns)
        self.assertIn("meta2", results["seq_len_df"].columns)
        self.assertEqual(results["metadata_columns"], ["meta1", "meta2"])
        self.assertEqual(
            results["seq_len_df"]
            .loc[results["seq_len_df"]["sample"] == "s1", "meta1"]
            .iloc[0],
            "val1",
        )

    def test_generate_plotting_data_empty_input_dir(self):
        contigs = ContigSequencesDirFmt()
        results = generate_plotting_data(contigs, metadata_df=None, n_cpus=1)
        self.assertTrue(results["seq_gc_df"].empty)
        self.assertTrue(results["seq_len_df"].empty)
        self.assertTrue(results["cumulative_df"].empty)
        self.assertTrue(results["nx_df"].empty)
        self.assertEqual(results["metadata_columns"], [])

    @mock.patch(
        "q2_assembly.contig_qc.qc.render_spec", return_value="mocked_spec_json_meta"
    )
    @mock.patch("q2_assembly.contig_qc.qc.q2templates.render")
    def test_evaluate_contigs_with_metadata(
        self,
        mock_q2render,
        mock_render_spec_func,
    ):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(f"{tmpdir}/output_viz_meta")
            os.makedirs(output_dir, exist_ok=True)
            # os.makedirs(f"{output_dir}/data", exist_ok=True)

            # Metadata
            metadata_df = pd.DataFrame(
                {"colA": ["group1", "group2"], "colB": [100, 200]},
                index=pd.Index(["s1_eval", "s2_eval"], name="id"),
            )
            qiime2_metadata = Metadata(metadata_df)

            evaluate_contigs(
                str(output_dir), contigs, metadata=qiime2_metadata, n_cpus=1
            )

            # Assertions
            mock_q2render.assert_called_once()

            context = mock_q2render.call_args[1]["context"]

            # Check Vega specs (assert they are non-empty strings)
            self.assertTrue(
                isinstance(context["vega_contig_length_spec"], str)
                and len(context["vega_contig_length_spec"]) > 0
            )
            self.assertTrue(
                isinstance(context["vega_gc_content_spec"], str)
                and len(context["vega_gc_content_spec"]) > 0
            )
            self.assertTrue(
                isinstance(context["vega_cumulative_length_spec"], str)
                and len(context["vega_cumulative_length_spec"]) > 0
            )
            self.assertTrue(
                isinstance(context["vega_nx_curve_spec"], str)
                and len(context["vega_nx_curve_spec"]) > 0
            )

            # Check sample_metrics
            expected_sample_metrics = [
                {
                    "sample": "s1_eval",
                    "count": 2,
                    "mean": 40.0,
                    "longest": 50,
                    "n50": 50,
                    "n90": 30,
                    "l50": 1,
                    "l90": 2,
                    "total_length": 80.0,
                    "colA": "group1",
                    "colB": 100,
                },
                {
                    "sample": "s2_eval",
                    "count": 1,
                    "mean": 100.0,
                    "longest": 100,
                    "n50": 100,
                    "n90": 100,
                    "l50": 1,
                    "l90": 1,
                    "total_length": 100.0,
                    "colA": "group2",
                    "colB": 200,
                },
            ]
            observed_sample_metrics = sorted(
                json.loads(context["sample_metrics"]), key=lambda x: x["sample"]
            )
            self.assertEqual(observed_sample_metrics, expected_sample_metrics)

            # Check categories
            self.assertEqual(json.loads(context["categories"]), ["colA", "colB"])

            # Check metadata values
            expected_values_json = {
                "colA": ["group1", "group2"],
                "colB": ["100.0", "200.0"],
            }
            self.assertDictEqual(json.loads(context["values"]), expected_values_json)

            # Check sample_ids_by_metadata
            expected_sample_ids_json = {
                "all_samples": sorted(["s1_eval", "s2_eval"]),
                "colA": {"group1": ["s1_eval"], "group2": ["s2_eval"]},
                "colB": {"100.0": ["s1_eval"], "200.0": ["s2_eval"]},
            }
            observed_sample_ids_by_metadata = json.loads(
                context["sample_ids_by_metadata"]
            )
            # Sort the 'all_samples' list for consistent comparison
            if "all_samples" in observed_sample_ids_by_metadata:
                observed_sample_ids_by_metadata["all_samples"] = sorted(
                    observed_sample_ids_by_metadata["all_samples"]
                )
            self.assertDictEqual(
                observed_sample_ids_by_metadata, expected_sample_ids_json
            )

            # Check if data files were created (optional, but good for sanity)
            data_path = output_dir / "data"
            self.assertTrue((data_path / "sample_metrics.tsv").exists())
            self.assertTrue((data_path / "gc_content_data.arrow").exists())
            self.assertTrue((data_path / "contig_length_data.arrow").exists())
            self.assertTrue((data_path / "cumulative_length_data.arrow").exists())
            self.assertTrue((data_path / "nx_curve_data.arrow").exists())
