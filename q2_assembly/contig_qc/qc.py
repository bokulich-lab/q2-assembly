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
MAX_CUMULATIVE_POINTS_DEFAULT = 500


def _process_single_fasta(fasta_file_path: Path):
    """
    Processes a single FASTA file to extract contig statistics.

    Args:
        fasta_file_path (Path): Path to the input FASTA file.

    Returns:
        tuple: A tuple containing the sample ID (derived from the filename)
               and a dictionary with raw contig data:
               - "lengths": List of contig lengths.
               - "total_length": Sum of all contig lengths.
               - "num_contigs": Total number of contigs.
               - "gc": List of GC content percentages for each contig.
               - "sorted_lengths": List of contig lengths, sorted in descending order.
    """
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


def _get_gc_content(sample_id: str, raw_data: dict) -> list:
    """
    Prepares GC content data rows for a single sample.

    Args:
        sample_id (str): The identifier for the sample.
        raw_data (dict): Raw data for the sample. Expected keys: "lengths", "gc".

    Returns:
        list: List of dicts, each with "sample" and "gc" keys.
    """
    return [
        {"sample": sample_id, "gc": gc}
        for _, gc in zip(raw_data["lengths"], raw_data["gc"])
    ]


def _get_contig_lengths(sample_id: str, raw_data: dict) -> list:
    """
    Prepares contig length data rows for a single sample.

    Args:
        sample_id (str): The identifier for the sample.
        raw_data (dict): Raw data for the sample. Expected keys: "lengths", "gc".

    Returns:
        list: List of dicts, each with "sample" and "contig_length" keys.
    """
    return [
        {"sample": sample_id, "contig_length": l}
        for l, _ in zip(raw_data["lengths"], raw_data["gc"])
    ]


def _calculate_cumulative_length(
    sample_id: str,
    sorted_lengths: list,
    max_points: int = MAX_CUMULATIVE_POINTS_DEFAULT,
) -> pd.DataFrame | None:
    """
    Calculates a (potentially pruned) cumulative length DataFrame for a sample.

    If the number of contigs exceeds `max_points`, the data points are
    subsampled to `max_points` for performance.

    Args:
        sample_id (str): The identifier for the sample.
        sorted_lengths (list): List of contig lengths, sorted in descending order.
        max_points (int): The maximum number of data points to include in the
                          output DataFrame. Defaults to MAX_CUMULATIVE_POINTS_DEFAULT.

    Returns:
        pd.DataFrame | None: A DataFrame with "sample", "rank", and
                              "cumulative_length" columns, or None if no sorted_lengths.
    """
    if not sorted_lengths:
        return None

    cum_sum_full = np.cumsum(sorted_lengths)
    ranks_full = np.arange(1, len(sorted_lengths) + 1)  # Ranks are 1-based
    num_contigs = len(sorted_lengths)

    if num_contigs <= max_points:
        selected_ranks = ranks_full
        selected_cum_sum = cum_sum_full
    else:
        indices_float = np.linspace(0, num_contigs - 1, num=max_points)
        selected_indices = np.unique(indices_float.astype(int))
        selected_ranks = ranks_full[selected_indices]
        selected_cum_sum = cum_sum_full[selected_indices]

    return pd.DataFrame(
        {
            "sample": sample_id,
            "rank": selected_ranks,
            "cumulative_length": selected_cum_sum,
        }
    )


def _calculate_nx_metrics(
    sample_id: str, sorted_lengths: list, total_length: int
) -> list:
    """
    Calculates N(x) values (e.g., N50, N90) for a single sample.

    Args:
        sample_id (str): The identifier for the sample.
        sorted_lengths (list): List of contig lengths, sorted in descending order.
        total_length (int): The total summed length of all contigs for the sample.

    Returns:
        list: A list of dictionaries, each representing an N(x) value with
              "sample", "percent", and "nx" keys. Returns an empty list if
              total_length is 0 or no sorted_lengths.
    """
    nx_rows = []
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
    return nx_rows


def _calculate_all_metrics(sample_id: str, raw_data: dict):
    """
    Calculates all derived metrics for a single sample from its raw contig data.

    This includes data for GC content plots, contig length plots,
    cumulative length plots (with pruning), and N(x) curves.

    Args:
        sample_id (str): The identifier for the sample.
        raw_data (dict): A dictionary containing raw data for the sample,
                         as generated by `_process_single_fasta`.
                         Expected keys: "lengths", "gc", "sorted_lengths", "total_length".

    Returns:
        dict: A dictionary containing various data structures for plotting:
              - "seq_gc_rows": List of dicts for GC content per contig.
              - "seq_len_rows": List of dicts for contig lengths.
              - "cumulative_df_part": DataFrame for cumulative length plot
                                      (points may be pruned).
              - "nx_rows": List of dicts for N(x) values.
    """
    seq_gc_rows = _get_gc_content(sample_id, raw_data)
    seq_len_rows = _get_contig_lengths(sample_id, raw_data)

    sorted_lengths = raw_data["sorted_lengths"]
    total_length = raw_data["total_length"]

    cumulative_df_part = _calculate_cumulative_length(sample_id, sorted_lengths)
    # The MAX_CUMULATIVE_POINTS constant is now a default in the helper function

    nx_rows = _calculate_nx_metrics(sample_id, sorted_lengths, total_length)

    return {
        "seq_gc_rows": seq_gc_rows,
        "seq_len_rows": seq_len_rows,
        "cumulative_df_part": cumulative_df_part,
        "nx_rows": nx_rows,
    }


def generate_plotting_data(
    contigs_dir: str, metadata_df: pd.DataFrame = None, n_cpus: int = 1
):
    """
    Generates all data required for plotting the QC results.

    This function processes all FASTA files in the `contigs_dir`,
    calculates raw and derived metrics in parallel, and merges
    with metadata if provided.

    Args:
        contigs_dir (str): Path to the directory containing contig FASTA files
                           (e.g., ContigSequencesDirFmt).
        metadata_df (pd.DataFrame, optional): DataFrame containing sample metadata.
                                             Defaults to None.
        n_cpus (int, optional): Number of CPUs to use for parallel processing.
                                Defaults to 1.

    Returns:
        dict: A dictionary containing DataFrames and lists for visualization:
              - "seq_gc_df": DataFrame with GC content per contig, merged with metadata.
              - "seq_len_df": DataFrame with contig lengths, merged with metadata.
              - "cumulative_df": DataFrame for cumulative length plots, merged with metadata.
              - "nx_df": DataFrame for N(x) plots, merged with metadata.
              - "metadata_columns": List of metadata column names that were merged.
    """
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
            _calculate_all_metrics, args_for_derived_metrics
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
    """
    Computes summary metrics for each sample based on contig length data.

    Metrics include count, mean length, longest contig, N50, N90, and total length.
    Metadata categories are included in the output if provided.

    Args:
        contigs_data_df (pd.DataFrame): DataFrame containing contig lengths
                                        (typically 'seq_len_df' from
                                        `generate_plotting_data`), potentially
                                        merged with metadata. Must contain 'sample'
                                        and 'contig_length' columns.
        categories (List[str]): A list of metadata category names to include
                                in the output metrics for each sample.

    Returns:
        list: A list of dictionaries, where each dictionary represents a sample
              and its computed metrics, including any specified metadata categories.
    """
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
    """
    Removes default Bootstrap 3 CSS and JS files from the q2template assets.

    This is a workaround to allow using Bootstrap 5 styles which are
    manually adjusted in the HTML templates.

    Args:
        output_dir (str): The root directory of the QIIME 2 visualization.
    """
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
    """
    Renders a Vega-Lite JSON specification from a Jinja2 template.

    Args:
        template_name (str): The filename of the Jinja2 template for the Vega-Lite spec
                             (e.g., "contig_length_spec.json.j2").
        **kwargs: Arbitrary keyword arguments to pass to the Jinja2 template rendering.

    Returns:
        str: A string containing the rendered Vega-Lite JSON specification.
    """
    # Render the four separate Vega-Lite specs for the dashboard
    spec_template_fp = os.path.join(TEMPLATES, "contig_qc", "vega", template_name)
    with open(spec_template_fp) as f:
        spec_template = jinja2.Template(f.read())
    return spec_template.render(**kwargs)


def estimate_column_count(samples: set[str]) -> int:
    """
    Estimates the optimal number of columns for multi-panel Vega plots
    based on the maximum length of sample IDs.

    Aims to prevent sample IDs from being truncated in plot titles.

    Args:
        samples (set[str]): A set of unique sample IDs.

    Returns:
        int: The estimated number of columns (2, 3, or 4).
    """
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
    """
    Generates a QIIME 2 visualization for contig quality control.

    This function orchestrates the processing of contig data, calculation of
    various metrics, generation of data for plots, and rendering of the
    HTML visualization.

    Args:
        output_dir (str): The directory where the visualization output will be saved.
        contigs (ContigSequencesDirFmt): QIIME 2 artifact containing per-sample
                                         contig sequences.
        metadata (Metadata, optional): Sample metadata to be incorporated into
                                       the plots and tables. Defaults to None.
        n_cpus (int, optional): Number of CPUs to use for parallel data processing.
                                Defaults to 1.
    """
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
