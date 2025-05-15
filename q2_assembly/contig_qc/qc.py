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
from typing import List

import jinja2
import numpy as np
import pandas as pd
import panel as pn
from qiime2 import Metadata
import pkg_resources
import q2templates
from q2_types.per_sample_sequences import ContigSequencesDirFmt
from skbio import DNA, read


TEMPLATES = pkg_resources.resource_filename("q2_assembly", "assets")


def generate_plotting_data(contigs_dir, metadata=None):
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

    # Build per-sequence DataFrames.
    seq_data_gc, seq_data_len = [], []
    for sample, d in samples_data.items():
        for L, gc in zip(d['lengths'], d['gc']):
            seq_data_gc.append({"sample": sample, "gc": gc})
            seq_data_len.append({"sample": sample, "contig_length": L})
    seq_gc_df = pd.DataFrame(seq_data_gc)
    seq_len_df = pd.DataFrame(seq_data_len)

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

    if metadata is not None:
        metadata = metadata.reset_index().rename(columns={metadata.index.name or 'index': 'sample'})
        seq_gc_df = seq_gc_df.merge(metadata, on='sample', how='left')
        seq_len_df = seq_len_df.merge(metadata, on='sample', how='left')
        cumulative_df = cumulative_df.merge(metadata, on='sample', how='left')
        nx_df = nx_df.merge(metadata, on='sample', how='left')

    return {
        "seq_gc_df": seq_gc_df,
        "seq_len_df": seq_len_df,
        "cumulative_df": cumulative_df,
        "nx_df": nx_df,
        "metadata_columns": list(metadata.columns.drop('sample')) if metadata is not None else []
    }


def compute_sample_metrics(seq_gc_df: pd.DataFrame, categories: List[str]):
    metrics = []
    for sample, group in seq_gc_df.groupby("sample"):
        contig_lengths = group['contig_length'].tolist()
        contig_lengths.sort(reverse=True)
        total_length = sum(contig_lengths)
        count = len(contig_lengths)
        n50 = n90 = 0
        acc = 0
        for l in contig_lengths:
            acc += l
            if not n50 and acc >= total_length * 0.5:
                n50 = l
            if not n90 and acc >= total_length * 0.9:
                n90 = l
        metrics.append({
            'sample': sample,
            'count': count,
            'n50': n50,
            'n90': n90,
            'total_length': total_length,
            **{category: group[category].iloc[0] for category in categories}
        })
    return metrics


def _cleanup_bootstrap(output_dir):
    # Remove unwanted files
    # until Bootstrap 3 is replaced with v5, remove the v3 scripts as
    # the HTML files are adjusted to work with v5
    os.remove(
        os.path.join(
            output_dir, "q2templateassets", "css", "bootstrap.min.css"
        )
    )
    os.remove(
        os.path.join(
            output_dir, "q2templateassets", "js", "bootstrap.min.js"
        )
    )

def evaluate_contigs(
        output_dir: str,
        contigs: ContigSequencesDirFmt,
        metadata: Metadata = None
):
    metadata_df = metadata.to_dataframe()
    data = generate_plotting_data(contigs, metadata.to_dataframe())

    categories = metadata_df.columns.tolist()
    values = {
        x: metadata_df[x].dropna().unique().tolist() for x in metadata_df.columns
    }

    sample_metrics = compute_sample_metrics(data['seq_len_df'], categories)

    # Render the four separate Vega-Lite specs for the dashboard
    def render_spec(template_name, **kwargs):
        spec_template_fp = os.path.join(TEMPLATES, "contig_qc", template_name)
        with open(spec_template_fp) as f:
            spec_template = jinja2.Template(f.read())
        return spec_template.render(**kwargs)

    vega_contig_length_spec = render_spec(
        "vega_contig_length_spec.json.j2",
        seq_df=json.dumps(data['seq_len_df'].to_dict(orient='records')),
    )
    vega_nx_curve_spec = render_spec(
        "vega_nx_curve_spec.json.j2",
        nx_df=json.dumps(data['nx_df'].to_dict(orient='records')),
    )
    vega_gc_content_spec = render_spec(
        "vega_gc_content_spec.json.j2",
        seq_df=json.dumps(data['seq_gc_df'].to_dict(orient='records')),
    )
    vega_cumulative_length_spec = render_spec(
        "vega_cumulative_length_spec.json.j2",
        cumulative_df=json.dumps(data['cumulative_df'].to_dict(orient='records')),
    )

    templates = [
        os.path.join(TEMPLATES, "contig_qc", "index.html"),
        os.path.join(TEMPLATES, "contig_qc", "grouped.html")
    ]
    context = {
        "tabs": [
            {"title": "Sample metrics", "url": "index.html"},
            {"title": "Group metrics", "url": "grouped.html"}
        ],
        "vega_contig_length_spec": vega_contig_length_spec,
        "vega_nx_curve_spec": vega_nx_curve_spec,
        "vega_gc_content_spec": vega_gc_content_spec,
        "vega_cumulative_length_spec": vega_cumulative_length_spec,
        "sample_metrics": json.dumps(sample_metrics),
        "categories": json.dumps(categories),
        "values": json.dumps(values),
    }

    for d in ("js", "css"):
        shutil.copytree(
            os.path.join(TEMPLATES, "contig_qc", d),
            os.path.join(output_dir, d)
        )

    q2templates.render(templates, output_dir, context=context)

    _cleanup_bootstrap(output_dir)