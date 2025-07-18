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
from importlib import resources
from multiprocessing import Pool
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.ipc as ipc
import q2templates
import qiime2
from q2_types.per_sample_sequences import ContigSequencesDirFmt
from qiime2 import Metadata, Artifact
from skbio import DNA, read

TEMPLATES = resources.files("q2_assembly") / "assets"
MAX_CUMULATIVE_POINTS = 500


def _process_single_fasta(fp: Path):
    """
    Processes a single FASTA file to extract contig statistics.

    Args:
        fp (Path): Path to the input FASTA file.

    Returns:
        tuple: A tuple containing the sample ID (derived from the filename)
               and a dictionary with raw contig data:
               - "lengths": List of contig lengths.
               - "total_length": Sum of all contig lengths.
               - "num_contigs": Total number of contigs.
               - "gc": List of GC content percentages for each contig.
               - "sorted_lengths": List of contig lengths, sorted in descending order.
    """
    sample_id = fp.stem.replace("_contigs", "")

    gc_vals, lengths = [], []

    for seq in read(str(fp), format="fasta", constructor=DNA):
        lengths.append(len(seq))
        gc_vals.append(100 * seq.gc_content() if len(seq) > 0 else 0.0)

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
    return [{"sample": sample_id, "gc": gc} for gc in raw_data["gc"]]


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
        {"sample": sample_id, "contig_length": contig_len}
        for contig_len in raw_data["lengths"]
    ]


def _calculate_cumulative_length(
    sample_id: str,
    sorted_lengths: list,
    max_points: int = MAX_CUMULATIVE_POINTS,
) -> pd.DataFrame | None:
    """
    Calculates a (potentially pruned) cumulative length DataFrame for a sample.

    If the number of contigs exceeds `max_points`, the data points are
    subsampled to `max_points` for performance.

    Args:
        sample_id (str): The identifier for the sample.
        sorted_lengths (list): List of contig lengths, sorted in descending order.
        max_points (int): The maximum number of data points to include in the
                          output DataFrame. Defaults to MAX_CUMULATIVE_POINTS.

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
    Calculates N(x) and L(x) values (e.g., N50, L50) for a single sample.

    Args:
        sample_id (str): The identifier for the sample.
        sorted_lengths (list): List of contig lengths, sorted in descending order.
        total_length (int): The total summed length of all contigs for the sample.

    Returns:
        list: A list of dictionaries, each representing an N(x)/L(x) value with
              "sample", "percent", "nx", and "lx" keys. Returns an empty list if
              total_length is 0 or no sorted_lengths.
    """
    nx_rows = []
    if total_length > 0 and sorted_lengths:
        cum_sum_nx = np.cumsum(sorted_lengths)
        percents = np.arange(1, 101)
        thresholds = total_length * percents / 100.0
        # Find the indices of the contigs where cumulative sum meets the threshold
        # idxs[i] is the 0-based index into sorted_lengths for percents[i]
        idxs = np.searchsorted(cum_sum_nx, thresholds, side="left").clip(
            max=len(sorted_lengths) - 1
        )
        nx_rows = [
            {
                "sample": sample_id,
                "percent": int(p),
                "nx": int(sorted_lengths[i]),
                "lx": i + 1,
            }
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
                         Expected keys: "lengths", "gc", "sorted_lengths",
                         "total_length".

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
    nx_rows = _calculate_nx_metrics(sample_id, sorted_lengths, total_length)

    return {
        "seq_gc_rows": seq_gc_rows,
        "seq_len_rows": seq_len_rows,
        "cumulative_df_part": cumulative_df_part,
        "nx_rows": nx_rows,
    }


def _df_to_arrow(df: pd.DataFrame, filename: str, output_dir: str):
    """
    Writes a pandas DataFrame to an Arrow IPC file.

    Args:
        df (pd.DataFrame): The DataFrame to write.
        filename (str): The name of the Arrow file to create.
        output_dir (str): The directory where the 'data' subdirectory
                          will be located or created, and the Arrow file
                          will be saved.
    """
    table = pa.Table.from_pandas(df)
    fp = os.path.join(output_dir, "data", filename)
    with ipc.RecordBatchFileWriter(fp, table.schema) as writer:
        writer.write_table(table)


def _reset_indices(artifacts: List[Artifact]) -> Artifact:
    """
    Concatenates a list of QIIME 2 Metadata Artifacts, resets their indices,
    and returns a new Artifact.

    Args:
        artifacts (List[Artifact]): A list of QIIME 2 Artifacts, each viewable as
                             qiime2.Metadata containing a pandas DataFrame.

    Returns:
        Artifact: A new QIIME 2 Artifact of type 'ImmutableMetadata' containing
                  the concatenated and re-indexed DataFrame.
    """
    _df = pd.concat([x.view(qiime2.Metadata).to_dataframe() for x in artifacts])
    _df.reset_index(inplace=True, drop=True)
    _df.index = _df.index.map(str)
    _df.index.name = "id"
    return qiime2.Artifact.import_data("ImmutableMetadata", qiime2.Metadata(_df))


def dump_all_to_arrow(data: dict, output_dir: str):
    """
    Dumps multiple DataFrames from a dictionary to Arrow IPC files.

    Creates a 'data' subdirectory in `output_dir` if it doesn't exist.
    The DataFrames are expected under keys 'seq_len_df', 'seq_gc_df',
    'cumulative_df', and 'nx_df' in the `data` dictionary.

    Args:
        data (dict): A dictionary containing pandas DataFrames.
        output_dir (str): The base directory where the 'data' subdirectory
                          will be created and Arrow files stored.
    """
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
        _df_to_arrow(df, fn, output_dir)


def generate_plotting_data(
    contigs_dir: ContigSequencesDirFmt,
    n_cpus: int = 1,
):
    """
    Generates all data required for plotting the QC results.

    This function processes all FASTA files in the `contigs_dir`,
    calculates raw and derived metrics in parallel.

    Args:
        contigs_dir (ContigSequencesDirFmt):Directory containing contig FASTA files.
        n_cpus (int, optional): Number of CPUs to use for parallel processing.
                                Defaults to 1.

    Returns:
        dict: A dictionary containing DataFrames for visualization:
              - "seq_gc_df": DataFrame with GC content per contig.
              - "seq_len_df": DataFrame with contig lengths per contig.
              - "cumulative_df": DataFrame for cumulative length plots.
              - "nx_df": DataFrame for N(x) plots.
    """
    fasta_files = list(Path(str(contigs_dir)).glob("*.fa"))

    all_seq_gc_rows = []
    all_seq_len_rows = []
    all_cumulative_df_parts = []
    all_nx_rows = []

    with Pool(processes=n_cpus) as pool:
        # Step 1: Process individual FASTA files to get raw data
        raw_sample_results = pool.map(_process_single_fasta, fasta_files)

        # Step 2: Calculate all derived metrics per sample in parallel
        args_for_derived_metrics = [
            (sample_id, data_dict)
            for sample_id, data_dict in raw_sample_results
            if data_dict
        ]
        derived_metrics_results = pool.starmap(
            _calculate_all_metrics, args_for_derived_metrics
        )

        # Step 3: Aggregate results from the second parallel step
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
        else pd.DataFrame(columns=["sample", "percent", "nx", "lx"])
    )

    for df in (seq_gc_df, seq_len_df, cumulative_df, nx_df):
        df.index = df.index.map(str)
        df.index.name = "id"

    return {
        "seq_gc_df": seq_gc_df,
        "seq_len_df": seq_len_df,
        "cumulative_df": cumulative_df,
        "nx_df": nx_df,
    }


def compute_sample_metrics(contigs_data_df: pd.DataFrame) -> pd.DataFrame:
    """
    Computes summary metrics for each sample based on contig length data.

    Metrics include count, mean length, longest contig, N50, N90, L50, L90,
    and total length.

    Args:
        contigs_data_df (pd.DataFrame): DataFrame containing contig lengths
                                        (typically 'seq_len_df' from
                                        `generate_plotting_data`). Must contain 'sample'
                                        and 'contig_length' columns.

    Returns:
        pd.DataFrame: DataFrame with each row representing a sample and its computed
                      metrics. Index is 'id', columns include 'sample' and metrics.
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
        l50 = l90 = 0
        nx = lx = 0
        for contig_len in contig_lengths:
            nx += contig_len
            lx += 1
            if not n50 and nx >= total_length * 0.5:
                n50 = contig_len
                l50 = lx
            if not n90 and nx >= total_length * 0.9:
                n90 = contig_len
                l90 = lx
        metrics.append(
            {
                "sample": sample,
                "count": count,
                "mean": mean,
                "longest": longest,
                "n50": n50,
                "n90": n90,
                "l50": l50,
                "l90": l90,
                "total_length": total_length,
            }
        )
    result = pd.DataFrame(metrics)
    if result.empty:
        result = pd.DataFrame(
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
            ]
        ).set_index("sample", drop=False)
    else:
        result = result.set_index("sample", drop=False)
    result.sort_index(inplace=True)
    result.index.name = "id"
    return result


def process_metadata(
    metadata: qiime2.Metadata,
    samples_by_metadata: dict,
) -> (dict, dict):
    """
    Processes categorical metadata for use in the visualization.

    Filters for categorical columns, fills NA values with the string "NA",
    and structures the data for dropdown menus and sample filtering in
    the HTML visualization.

    Args:
        metadata (qiime2.Metadata): The QIIME 2 metadata object.
        samples_by_metadata (dict): A dictionary to be populated with sample IDs
                                    grouped by metadata category and value.
                                    It's typically initialized with an
                                    `all_samples` key.

    Returns:
        tuple:
            - dict: A dictionary where keys are metadata categories (column names)
                    and values are lists of unique stringified values for that
                    category (e.g., `{\"category1\": [\"valueA\", \"valueB\"]}`).
            - dict: The `samples_by_metadata` dictionary, updated to include
                    mappings from each category and its values to the list of
                    associated sample IDs (e.g., `{\"category1\": {\"valueA\":
                    [\"sample1\", \"sample2\"], ...}}`).
    """
    metadata = metadata.filter_columns(column_type="categorical").to_dataframe()
    metadata.fillna("NA")
    metadata_values = {x: metadata[x].unique().tolist() for x in metadata.columns}
    for cat, values in metadata_values.items():
        samples_by_metadata[cat] = {}
        for cat_value in values:
            cat_samples = metadata[metadata[cat] == cat_value].index.unique().tolist()
            samples_by_metadata[cat][cat_value] = sorted(cat_samples)
    return metadata_values, samples_by_metadata


def _evaluate_contigs(
    contigs: ContigSequencesDirFmt,
    n_cpus: int = 1,
) -> (
    qiime2.Metadata,
    qiime2.Metadata,
    qiime2.Metadata,
    qiime2.Metadata,
    qiime2.Metadata,
):
    data = generate_plotting_data(contigs, n_cpus=n_cpus)
    sample_metrics = compute_sample_metrics(data["seq_len_df"])
    return (
        qiime2.Metadata(sample_metrics),
        qiime2.Metadata(data["nx_df"]),
        qiime2.Metadata(data["seq_gc_df"]),
        qiime2.Metadata(data["seq_len_df"]),
        qiime2.Metadata(data["cumulative_df"]),
    )


def _visualize_contig_qc(
    output_dir: str,
    per_sample_metrics: qiime2.Metadata,
    nx: qiime2.Metadata,
    gc: qiime2.Metadata,
    lengths: qiime2.Metadata,
    cumulative: qiime2.Metadata,
    metadata: Metadata = None,
):
    data = {
        "seq_gc_df": gc.to_dataframe(),
        "seq_len_df": lengths.to_dataframe(),
        "cumulative_df": cumulative.to_dataframe(),
        "nx_df": nx.to_dataframe(),
        "per_sample": per_sample_metrics.to_dataframe(),
    }
    samples_by_metadata = {
        "all_samples": sorted(list(data["seq_len_df"]["sample"].unique()))
    }

    # Prepare metadata categories/values for the dropdowns
    # Merge metadata into all DataFrames for the Vega visualizations
    if metadata:
        metadata_context_values, samples_by_metadata = process_metadata(
            metadata, samples_by_metadata
        )
        metadata_df = (
            metadata.filter_columns(column_type="categorical")
            .to_dataframe()
            .fillna("NA")
        )
        for key, df in data.items():
            data[key] = df.merge(
                metadata_df, how="left", left_on="sample", right_index=True
            )
    else:
        metadata_context_values = {}

    # Dump all the data for Vega into arrow tables for more efficient in-browser access
    dump_all_to_arrow(data, output_dir)

    # Prepare templates and the context
    templates = [
        TEMPLATES / "contig_qc" / "index.html",
        TEMPLATES / "contig_qc" / "table.html",
    ]
    context = {
        "tabs": [
            {"title": "Sample metrics", "url": "index.html"},
            {"title": "Table view", "url": "table.html"},
        ],
        "has_metadata": "true" if metadata is not None else "false",
        "sample_metrics": json.dumps(data["per_sample"].to_dict("records")),
        "categories": json.dumps(list(metadata_context_values.keys())),
        "values": json.dumps(metadata_context_values),
        "sample_ids_by_metadata": json.dumps(samples_by_metadata),
        "page_size": 100,
    }
    context.update(
        {
            f"vega_{x}_spec": json.dumps(
                json.load(open(TEMPLATES / "contig_qc" / "vega" / f"{x}.json"))
            )
            for x in ("contig_length", "nx_curve", "gc_content", "cumulative_length")
        }
    )
    if metadata is not None:
        templates.insert(1, TEMPLATES / "contig_qc" / "grouped.html")
        context["tabs"].insert(1, {"title": "Group metrics", "url": "grouped.html"})

    # Dump per-sample metrics into a tsv file
    data["per_sample"].to_csv(
        os.path.join(output_dir, "data", "sample_metrics.tsv"), sep="\t", index=True
    )

    # Copy JS/CSS files
    for d in ("js", "css"):
        shutil.copytree(TEMPLATES / "contig_qc" / d, os.path.join(output_dir, d))

    q2templates.render(templates, output_dir, context=context)


def evaluate_contigs(ctx, contigs, metadata=None, n_cpus=1, num_partitions=1):
    _partition = ctx.get_action("assembly", "partition_contigs")
    _evaluate = ctx.get_action("assembly", "_evaluate_contigs")
    _visualize = ctx.get_action("assembly", "_visualize_contig_qc")

    # Do not partition if only 1 partition was requested
    if num_partitions == 1:
        (
            per_sample_results,
            nx,
            gc,
            lengths,
            cumulative,
        ) = _evaluate(contigs, n_cpus)
        (viz,) = _visualize(
            per_sample_results,
            nx,
            gc,
            lengths,
            cumulative,
            metadata,
        )
    else:
        (partitioned_contigs,) = _partition(contigs, num_partitions)
        results_all, nx_all, gc_all, lengths_all, cumulative_all = [], [], [], [], []
        for _contigs in partitioned_contigs.values():
            (
                per_sample,
                nx,
                gc,
                lengths,
                cumulative,
            ) = _evaluate(_contigs, n_cpus)

            results_all.append(per_sample)
            nx_all.append(nx)
            gc_all.append(gc)
            lengths_all.append(lengths)
            cumulative_all.append(cumulative)

        per_sample_results = ctx.make_artifact(
            "ImmutableMetadata",
            qiime2.Metadata(
                pd.concat([x.view(qiime2.Metadata).to_dataframe() for x in results_all])
            ),
        )

        artifacts = {
            "nx": nx_all,
            "gc": gc_all,
            "lengths": lengths_all,
            "cumulative": cumulative_all,
        }
        for key, _artifacts in artifacts.items():
            artifacts[key] = _reset_indices(_artifacts)

        (viz,) = _visualize(
            per_sample_results,
            artifacts["nx"],
            artifacts["gc"],
            artifacts["lengths"],
            artifacts["cumulative"],
            metadata,
        )

    return per_sample_results, viz
