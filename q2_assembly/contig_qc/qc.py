# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import os
import shutil
from multiprocessing import Pool
from pathlib import Path
from typing import List

import jinja2
import numpy as np
import pandas as pd
import pkg_resources
import pyarrow as pa
import pyarrow.ipc as ipc
import q2templates
from q2_types.per_sample_sequences import ContigSequencesDirFmt
from qiime2 import Metadata
from skbio import DNA, read

TEMPLATES = pkg_resources.resource_filename("q2_assembly", "assets")


def _process_single_fasta(fasta_file_path: Path):
    sample_id = fasta_file_path.stem.replace("_contigs", "")

    sequences = list(read(str(fasta_file_path), format="fasta", constructor=DNA))
    lengths = [len(seq) for seq in sequences]
    gc_vals = [100 * seq.gc_content() if len(seq) > 0 else 0.0 for seq in sequences]

    sample_data = {
        "lengths": lengths,
        "total_length": sum(lengths),
        "num_contigs": len(lengths),
        "gc": gc_vals,
        "sorted_lengths": sorted(lengths, reverse=True),
    }
    return sample_id, sample_data


def _calculate_all_derived_data_for_sample(sample_id: str, raw_data: dict):
    """Calculates all derived metrics for a single sample from its raw data."""
    seq_gc_rows = [
        {"sample": sample_id, "gc": gc}
        for _, gc in zip(raw_data["lengths"], raw_data["gc"])
    ]
    seq_len_rows = [
        {"sample": sample_id, "contig_length": l}
        for l, _ in zip(raw_data["lengths"], raw_data["gc"])
    ]

    cumulative_df_part = None
    sorted_lengths = raw_data["sorted_lengths"]
    # Max points for the cumulative length plot per sample.
    # Chosen to provide good visual detail without overwhelming the browser.
    MAX_CUMULATIVE_POINTS = 500

    if sorted_lengths:
        cum_sum_full = np.cumsum(sorted_lengths)
        # Ranks are 1-based for plotting
        ranks_full = np.arange(1, len(sorted_lengths) + 1)
        num_contigs = len(sorted_lengths)

        if num_contigs <= MAX_CUMULATIVE_POINTS:
            # If fewer or an equal number of contigs than max_points, use all data
            selected_ranks = ranks_full
            selected_cum_sum = cum_sum_full
        else:
            # Subsample the data to approximately MAX_CUMULATIVE_POINTS
            # np.linspace will ensure the first (0) and last (num_contigs-1) indices are covered.
            indices_float = np.linspace(
                0, num_contigs - 1, num=MAX_CUMULATIVE_POINTS
            )
            # Convert to integer indices and take unique ones (handles cases where
            # MAX_CUMULATIVE_POINTS is close to num_contigs, preventing duplicates
            # and ensuring sorted order).
            selected_indices = np.unique(indices_float.astype(int))

            selected_ranks = ranks_full[selected_indices]
            selected_cum_sum = cum_sum_full[selected_indices]

        cumulative_df_part = pd.DataFrame(
            {"sample": sample_id, "rank": selected_ranks, "cumulative_length": selected_cum_sum}
        )

    nx_rows = []
    total_length = raw_data["total_length"]
    if total_length > 0 and sorted_lengths:
        cum_sum_nx = np.cumsum(sorted_lengths)
        percents = np.arange(1, 101)
        thresholds = total_length * percents / 100.0
        idxs = np.searchsorted(cum_sum_nx, thresholds, side="left").clip(
            max=len(sorted_lengths) - 1
        )
        nx_rows = [
            {"sample": sample_id, "percent": int(p), "nx": int(sorted_lengths[i])}
            for p, i in zip(percents, idxs)
        ]

    return {
        "seq_gc_rows": seq_gc_rows,
        "seq_len_rows": seq_len_rows,
        "cumulative_df_part": cumulative_df_part,
        "nx_rows": nx_rows,
    }


def generate_plotting_data(
    contigs_dir: str, metadata_df: pd.DataFrame = None, n_cpus: int = 1
):
    fasta_files = list(Path(str(contigs_dir)).glob("*.fa"))

    all_seq_gc_rows = []
    all_seq_len_rows = []
    all_cumulative_df_parts = []
    all_nx_rows = []

    with Pool(processes=n_cpus) as pool:
        # Step 1: Process individual FASTA files to get raw data
        raw_sample_results = pool.map(_process_single_fasta, fasta_files)

        # Prepare arguments for the second parallel step
        args_for_derived_metrics = [
            (sample_id, data_dict)
            for sample_id, data_dict in raw_sample_results
            if data_dict
        ]

        # Step 2: Calculate all derived metrics per sample in parallel
        derived_metrics_results = pool.starmap(
            _calculate_all_derived_data_for_sample, args_for_derived_metrics
        )

        # Step 3: Aggregate results from the second parallel step (serial, but fast)
        for res_dict in derived_metrics_results:
            all_seq_gc_rows.extend(res_dict["seq_gc_rows"])
            all_seq_len_rows.extend(res_dict["seq_len_rows"])
            all_nx_rows.extend(res_dict["nx_rows"])
            if res_dict["cumulative_df_part"] is not None:
                all_cumulative_df_parts.append(res_dict["cumulative_df_part"])

    # Create final DataFrames
    seq_gc_df = (
        pd.DataFrame(all_seq_gc_rows)
        if all_seq_gc_rows
        else pd.DataFrame(columns=["sample", "gc"])
    )
    seq_len_df = (
        pd.DataFrame(all_seq_len_rows)
        if all_seq_len_rows
        else pd.DataFrame(columns=["sample", "contig_length"])
    )

    if all_cumulative_df_parts:
        cumulative_df = pd.concat(all_cumulative_df_parts, ignore_index=True)
    else:
        cumulative_df = pd.DataFrame(columns=["sample", "rank", "cumulative_length"])

    nx_df = (
        pd.DataFrame(all_nx_rows)
        if all_nx_rows
        else pd.DataFrame(columns=["sample", "percent", "nx"])
    )

    metadata_columns_list = []
    if metadata_df is not None:
        metadata_df = metadata_df.fillna("NA")

        seq_gc_df = seq_gc_df.merge(
            metadata_df, left_on="sample", right_index=True, how="left"
        )
        seq_len_df = seq_len_df.merge(
            metadata_df, left_on="sample", right_index=True, how="left"
        )
        cumulative_df = cumulative_df.merge(
            metadata_df, left_on="sample", right_index=True, how="left"
        )
        nx_df = nx_df.merge(metadata_df, left_on="sample", right_index=True, how="left")
        metadata_columns_list = list(
            metadata_df.columns.drop("sample", errors="ignore")
        )

    return {
        "seq_gc_df": seq_gc_df,
        "seq_len_df": seq_len_df,
        "cumulative_df": cumulative_df,
        "nx_df": nx_df,
        "metadata_columns": metadata_columns_list,
    }


def compute_sample_metrics(contigs_data_df: pd.DataFrame, categories: List[str]):
    metrics = []
    for sample, group in contigs_data_df.groupby("sample"):
        contig_lengths = group["contig_length"].tolist()
        contig_lengths.sort(reverse=True)
        total_length = sum(contig_lengths)
        count = len(contig_lengths)
        mean = np.round(np.mean(contig_lengths), 0)
        longest = contig_lengths[0]
        n50 = n90 = 0
        acc = 0
        for l in contig_lengths:
            acc += l
            if not n50 and acc >= total_length * 0.5:
                n50 = l
            if not n90 and acc >= total_length * 0.9:
                n90 = l
        metrics.append(
            {
                "sample": sample,
                "count": count,
                "mean": mean,
                "longest": longest,
                "n50": n50,
                "n90": n90,
                "total_length": total_length,
                **{category: group[category].iloc[0] for category in categories},
            }
        )
    return metrics


def _cleanup_bootstrap(output_dir):
    # Remove unwanted files
    # until Bootstrap 3 is replaced with v5, remove the v3 scripts as
    # the HTML files are adjusted to work with v5
    bootstrap_css_path = (
        Path(output_dir) / "q2templateassets" / "css" / "bootstrap.min.css"
    )
    bootstrap_js_path = (
        Path(output_dir) / "q2templateassets" / "js" / "bootstrap.min.js"
    )

    for fp in (bootstrap_css_path, bootstrap_js_path):
        try:
            os.remove(fp)
        except FileNotFoundError:
            pass  # File already removed or never existed


def render_spec(template_name, **kwargs):
    # Render the four separate Vega-Lite specs for the dashboard
    spec_template_fp = os.path.join(TEMPLATES, "contig_qc", "vega", template_name)
    with open(spec_template_fp) as f:
        spec_template = jinja2.Template(f.read())
    return spec_template.render(**kwargs)


def estimate_column_count(samples: set[str]) -> int:
    max_len = max([len(x) for x in samples])
    if max_len >= 16:
        return 2
    elif max_len >= 9:
        return 3
    else:
        return 4


def evaluate_contigs(
    output_dir: str,
    contigs: ContigSequencesDirFmt,
    metadata: Metadata = None,
    n_cpus: int = 1,
):
    metadata_df = None
    # Categories and values for the visualization context (dropdowns, table headers)
    # These come directly from the input metadata, if provided.
    metadata_context_categories = []
    metadata_context_values = {}

    if metadata:
        metadata_df = metadata.to_dataframe()
        metadata_context_categories = metadata_df.columns.tolist()
        metadata_context_values = {
            x: [str(y) for y in metadata_df[x].dropna().unique().tolist()]
            for x in metadata_df.columns
        }

    data = generate_plotting_data(contigs, metadata_df=metadata_df, n_cpus=n_cpus)

    # Categories for sample_metrics should be those actually merged into the dataframes
    # and available for grouping/display in the metrics table.
    # These are returned by generate_plotting_data as 'metadata_columns'.
    categories_for_metrics = data["metadata_columns"]
    sample_metrics = compute_sample_metrics(data["seq_len_df"], categories_for_metrics)

    n_cols = estimate_column_count(set(data["seq_len_df"]["sample"]))

    os.makedirs(os.path.join(output_dir, "data"))

    data_for_export = zip(
        [data["seq_len_df"], data["seq_gc_df"], data["cumulative_df"], data["nx_df"]],
        [
            "contig_length_data.arrow",
            "gc_content_data.arrow",
            "cumulative_length_data.arrow",
            "nx_curve_data.arrow",
        ],
    )
    for df, fn in data_for_export:
        table = pa.Table.from_pandas(df)
        fp = os.path.join(output_dir, "data", fn)
        with ipc.RecordBatchFileWriter(fp, table.schema) as writer:
            writer.write_table(table)

    # Prepare sample_ids_by_metadata for the custom legend
    sample_ids_by_metadata = {
        "all_samples": sorted(list(data["seq_len_df"]["sample"].unique()))
    }

    if metadata_df is not None:  # Check if metadata was provided
        for category_name in data["metadata_columns"]:
            sample_ids_by_metadata[category_name] = {}
            unique_values_for_category = (
                data["seq_len_df"][category_name].dropna().unique()
            )
            for cat_value in unique_values_for_category:
                relevant_samples = (
                    data["seq_len_df"][data["seq_len_df"][category_name] == cat_value][
                        "sample"
                    ]
                    .unique()
                    .tolist()
                )
                sample_ids_by_metadata[category_name][str(cat_value)] = sorted(
                    relevant_samples
                )

    templates = [
        os.path.join(TEMPLATES, "contig_qc", "index.html"),
        os.path.join(TEMPLATES, "contig_qc", "grouped.html"),
        os.path.join(TEMPLATES, "contig_qc", "table.html"),
    ]
    context = {
        "tabs": [
            {"title": "Sample metrics", "url": "index.html"},
            {"title": "Group metrics", "url": "grouped.html"},
            {"title": "Table view", "url": "table.html"},
        ],
        "vega_contig_length_spec": render_spec(
            "contig_length_spec.json.j2",
            n_cols=n_cols,
        ),
        "vega_nx_curve_spec": render_spec(
            "nx_curve_spec.json.j2",
            n_cols=n_cols,
        ),
        "vega_gc_content_spec": render_spec(
            "gc_content_spec.json.j2",
            n_cols=n_cols,
        ),
        "vega_cumulative_length_spec": render_spec(
            "cumulative_length_spec.json.j2",
            n_cols=n_cols,
        ),
        "sample_metrics": json.dumps(sample_metrics),
        "categories": json.dumps(metadata_context_categories),
        "values": json.dumps(metadata_context_values),
        "sample_ids_by_metadata": json.dumps(sample_ids_by_metadata),
        "page_size": 100,
    }

    pd.DataFrame(sample_metrics)[
        ["sample", "count", "mean", "n50", "n90", "total_length"]
    ].to_csv(os.path.join(output_dir, "sample_metrics.tsv"), sep="\t", index=False)

    for d in ("js", "css"):
        shutil.copytree(
            os.path.join(TEMPLATES, "contig_qc", d), os.path.join(output_dir, d)
        )

    q2templates.render(templates, output_dir, context=context)

    _cleanup_bootstrap(output_dir)
