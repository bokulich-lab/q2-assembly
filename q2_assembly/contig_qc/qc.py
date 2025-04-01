# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import os
import glob
import re
import shutil

import numpy as np
import pandas as pd
import pkg_resources
import q2templates
from q2_types.per_sample_sequences import ContigSequencesDirFmt
from skbio import DNA, read
import altair as alt


TEMPLATES = pkg_resources.resource_filename("q2_assembly", "assets")


def generate_plotting_data(contigs_dir):
    # Process each FASTA file and build sample-level info.
    samples_data = {}
    for file in glob.glob(os.path.join(str(contigs_dir), "*.fa")):
        basename = os.path.basename(file)
        # Extract sample name from pattern: "<sample id>_contigs.fa"
        m = re.match(r"(.*)_contigs\.fa$", basename)
        sample = m.group(1) if m else basename
        sequences = list(read(file, format="fasta", constructor=DNA))
        lengths = [len(seq) for seq in sequences]
        gc_vals = [100 * (str(seq).count('G') + str(seq).count('C')) / (len(seq) if len(seq) > 0 else 1)
                   for seq in sequences]
        samples_data[sample] = {
            'lengths': lengths,
            'total_length': sum(lengths),
            'num_contigs': len(lengths),
            'gc': gc_vals,
            'sorted_lengths': sorted(lengths, reverse=True)
        }

    # Build a per-sequence DataFrame.
    seq_data = []
    for sample, d in samples_data.items():
        for L, gc in zip(d['lengths'], d['gc']):
            seq_data.append({"sample": sample, "contig_length": L, "gc": gc})
    seq_df = pd.DataFrame(seq_data)

    # Compute cumulative contig lengths per sample.
    cum_list = []
    for sample, d in samples_data.items():
        sl = d['sorted_lengths']
        if sl:
            cum = np.cumsum(sl)
            ranks = np.arange(1, len(cum) + 1)
            cum_list.append(pd.DataFrame({
                "sample": sample,
                "rank": ranks,
                "cumulative_length": cum
            }))
    cumulative_df = pd.concat(cum_list, ignore_index=True)

    # Compute Nx curves per sample.
    nx_rows = []
    for sample, d in samples_data.items():
        sl = d['sorted_lengths']
        total = d['total_length']
        if total > 0 and sl:
            cum = np.cumsum(sl)
            for p in range(1, 101):
                threshold = total * p / 100.0
                idx = np.searchsorted(cum, threshold, side="left")
                nx_val = sl[idx] if idx < len(sl) else sl[-1]
                nx_rows.append({"sample": sample, "percent": p, "nx": nx_val})
    nx_df = pd.DataFrame(nx_rows)

    return {
        "seq_df": seq_df,
        "cumulative_df": cumulative_df,
        "nx_df": nx_df
    }


def generate_visualization(data):
    # Define a multi-select for filtering via the legend.
    sample_sel = alt.selection_multi(fields=["sample"], bind="legend")
    # Define a click selection: single-click to highlight a sample.
    click_sel = alt.selection_single(fields=["sample"], on="click", empty="all", clear="dblclick")
    combined_sel = click_sel  # used for opacity highlighting

    # Define a dropdown selection parameter for sample filtering.
    sample_options = ["All"] + sorted(data["seq_df"]["sample"].unique().tolist())
    dropdown_sel = alt.param(
        name="dropdown",
        value="All",
        bind=alt.binding_select(options=sample_options, name="Select Sample:")
    )
    dropdown_filter = "dropdown === 'All' || datum.sample === dropdown"

    # 1. Contig Length Distribution as a Line Histogram.
    contig_hist_data = alt.Chart(data["seq_df"]).transform_bin(
        "binned", "contig_length", bin=alt.Bin(maxbins=30)
    ).transform_aggregate(
        count='count()', groupby=["binned", "sample"]
    )
    contig_hist = contig_hist_data.mark_line().encode(
        alt.X("binned:Q", axis=alt.Axis(format="~s", title="Contig Length")),
        alt.Y("count:Q", axis=alt.Axis(format="~s", title="Count")),
        alt.Color("sample:N", legend=alt.Legend()),
        alt.Tooltip("sample:N", title="Sample"),
        opacity=alt.condition(combined_sel, alt.value(1), alt.value(0.2))
    ).transform_filter(sample_sel
                       ).transform_filter(dropdown_filter
                                          ).add_params(dropdown_sel
                                                       ).add_selection(click_sel
                                                                       ).interactive().properties(
        width=500, height=300, title="Contig Length Distribution"
    )

    # 2. Nx Curve Plot.
    nx_curve = alt.Chart(data["nx_df"]).mark_line(interpolate="step-after").encode(
        alt.X("nx:Q", axis=alt.Axis(format="~s", title="Contig Length (Nx)")),
        alt.Y("percent:Q", title="Percent of Assembly"),
        alt.Color("sample:N", legend=alt.Legend()),
        alt.Tooltip(["sample:N", "nx:Q", "percent:Q"]),
        opacity=alt.condition(combined_sel, alt.value(1), alt.value(0.2))
    ).transform_filter(sample_sel
                       ).transform_filter(dropdown_filter
                                          ).add_params(dropdown_sel
                                                       ).add_selection(click_sel
                                                                       ).interactive().properties(
        width=500, height=300, title="Nx Curve per Sample"
    )

    # 3. GC Content Distribution as Density Lines.
    gc_violin = alt.Chart(data["seq_df"]).transform_density(
        "gc",
        as_=["gc", "density"],
        groupby=["sample"],
        extent=[0, 100]
    ).mark_line(strokeWidth=3).encode(
        alt.X("gc:Q", title="GC Content (%)"),
        alt.Y("density:Q", title="Density"),
        alt.Color("sample:N", legend=alt.Legend()),
        alt.Tooltip("sample:N", title="Sample"),
        opacity=alt.condition(combined_sel, alt.value(1), alt.value(0.2))
    ).transform_filter(sample_sel
                       ).transform_filter(dropdown_filter
                                          ).add_params(dropdown_sel
                                                       ).add_selection(click_sel
                                                                       ).interactive().properties(
        width=500, height=300, title="GC Content Distribution per Sample"
    )

    # 4. Cumulative Contig Length Curves.
    cumulative_chart = alt.Chart(data["cumulative_df"]).mark_line().encode(
        alt.X("rank:Q", title="No. of contigs"),
        alt.Y("cumulative_length:Q", axis=alt.Axis(format="~s", title="Cumulative Length")),
        alt.Color("sample:N", legend=alt.Legend()),
        alt.Tooltip("sample:N", title="Sample"),
        opacity=alt.condition(combined_sel, alt.value(1), alt.value(0.2))
    ).transform_filter(sample_sel
                       ).transform_filter(dropdown_filter
                                          ).add_params(dropdown_sel
                                                       ).add_selection(click_sel
                                                                       ).interactive().properties(
        width=500, height=300, title="Cumulative Contig Length Curves"
    )

    dashboard = alt.vconcat(
        alt.hconcat(contig_hist, nx_curve),
        alt.hconcat(gc_violin, cumulative_chart)
    ).configure_title(fontSize=16).add_selection(sample_sel)

    return dashboard


def evaluate_contigs_new(output_dir: str, contigs: ContigSequencesDirFmt):
    data = generate_plotting_data(contigs)
    dashboard = generate_visualization(data)

    shutil.copy(os.path.join(TEMPLATES, "contig_qc", "index.html"), output_dir)
    context = {
        'vega_json': json.dumps(dashboard.to_dict())
    }

    templates = [
        os.path.join(TEMPLATES, "contig_qc", "index.html"),
    ]

    q2templates.render(templates, output_dir, context=context)

