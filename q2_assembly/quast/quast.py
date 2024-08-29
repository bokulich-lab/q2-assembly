# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import json
import os
import platform
import subprocess
import tempfile
from distutils.dir_util import copy_tree
from typing import List, Union
from zipfile import ZipFile

import pandas as pd
import pkg_resources
import q2templates
from q2_types.feature_data import DNAFASTAFormat, DNAIterator
from q2_types.per_sample_sequences import (
    BAMDirFmt,
    ContigSequencesDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)

from q2_assembly.quast.utils import _parse_columns

from .._utils import (
    _construct_param,
    _get_sample_from_path,
    _modify_links,
    _process_common_input_params,
    _remove_html_element,
    run_command,
)

TEMPLATES = pkg_resources.resource_filename("q2_assembly", "assets")


def _process_quast_arg(arg_key, arg_val):
    """Creates a list with argument and its value.

    Argument values represented by a list will be converted to a single
    string joined by commas, e.g.: [1, 2, 3] -> '1,2,3'.
    Argument names will be converted to command line parameters by
    appending a '--' prefix and replacing all '_' with '-',
    e.g.: 'some_parameter' -> '--some-parameter'.

    Args:
        arg_key (str): Argument name.
        arg_val: Argument value.

    Returns:
        [converted_arg, arg_value]: List containing a prepared command line
            parameter and its value.
    """
    if isinstance(arg_val, bool) and arg_val:
        return [_construct_param(arg_key)]
    elif arg_key == "threads" and arg_val > 1 and platform.system() != "Linux":
        # This is a limitation in quast that is documented in this issue:
        # https://github.com/ablab/quast/issues/175
        raise ValueError("Multiprocessing is currently only supported on Linux.")
    elif not isinstance(arg_val, list):
        return [_construct_param(arg_key), str(arg_val)]
    else:
        arg_value = ",".join(str(x) for x in arg_val)
        return [_construct_param(arg_key), arg_value]


def _split_reference(ref: DNAFASTAFormat, all_refs_dir: str) -> List[str]:
    all_seq_fps = []
    for seq in ref.view(DNAIterator):
        seq_id = seq.metadata["id"]
        seq_fp = os.path.join(all_refs_dir, f"{seq_id}.fasta")
        all_seq_fps.append(seq_fp)
        with open(seq_fp, "w") as f:
            f.write(f">{seq_id}\n{str(seq)}")
    return all_seq_fps


def _evaluate_contigs(
    results_dir: str,
    contigs: ContigSequencesDirFmt,
    reads: dict,
    paired: bool,
    references: List[DNAFASTAFormat],
    mapped_reads: BAMDirFmt,
    common_args: list,
) -> List[str]:
    """Runs the contig assembly QC using QUAST.

    Constructs and runs the final QUAST command.

    Args:
        results_dir (str): Directory to store the final reports.
        contigs (ContigSequencesDirFmt): Contigs to be evaluated.
        reads (dict): Dictionary containing mapping of samples to their
            forward and reverse reads, e.g.:
            {'sample1': {'fwd': '/path/to/reads', 'rev': '/path/to/reads'}}.
        mapped_reads (BAMDirFmt): Mapping of reads to contigs
        common_args (list): List of common flags and their values for
            the QUAST command.

    Returns:
        result (list): List of sample names that were processed.
    """

    # TODO: this will probably get replaced by "quast.py" once we
    #  want to extend the support beyond metagenomes
    cmd = ["metaquast.py", "-o", results_dir]
    cmd.extend(common_args)
    samples = []

    if reads and mapped_reads:
        reads = None
        print("Both reads and mapped reads are provided. Reads will be ignored.")

    for fp in sorted(glob.glob(os.path.join(str(contigs), "*_contigs.fa"))):
        cmd.append(fp)
        samples.append(_get_sample_from_path(fp))

    if mapped_reads:
        bam_fps = sorted(glob.glob(os.path.join(str(mapped_reads), "*_alignment.bam")))
        cmd.extend(["--bam", ",".join(bam_fps)])
    elif reads:
        rev_count = sum([True if x["rev"] else False for _, x in reads.items()])
        # TODO: this is a strange statement which most likely
        #  doesn't do what it should - check and fix
        if (rev_count < len(samples) > rev_count) and paired:
            raise Exception(
                f"Number of reverse reads ({rev_count}) does not match "
                f"the number of provided contig files ({len(samples)}). "
                f"Please check your input files."
            )

        for s in samples:
            try:
                if paired:
                    cmd.extend(
                        ["--pe1", reads[s].get("fwd"), "--pe2", reads[s].get("rev")]
                    )
                else:
                    cmd.extend(["--single", reads[s].get("fwd")])
            except KeyError:
                # TODO: improve this msg by adding name of the missing sample
                raise Exception(
                    "Some samples are missing from the reads file. "
                    "Please check your input files."
                )

    if references:
        all_refs_dir = os.path.join(results_dir, "references")
        os.makedirs(all_refs_dir, exist_ok=True)
        all_ref_fps = []
        # we need to split the references into separate files so that QUAST
        # can correctly display alignment details per reference (otherwise it
        # will show those as if all the provided sequences belonged to a single
        # reference
        for ref in references:
            all_ref_fps.extend(_split_reference(ref, all_refs_dir))
        for fp in all_ref_fps:
            cmd.extend(["-r", fp])

    try:
        run_command(cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running "
            f"QUAST, (return code {e.returncode}), "
            "please inspect stdout and stderr to learn more."
        )

    return samples


def _fix_html_reports(results_dp: str):
    """Cleans up reports generated by QUAST.

    A few HTML elements generated by QUAST could break the visualization
    when clicked so they need to be removed or disabled.
    Modifies the files in-place.

    Args:
         results_dp (str): Directory where the reports are stored.

    """

    # remove link to icarus browser from report.html
    # (otherwise, the visualization can get messed up once clicked)
    report_fp = os.path.join(results_dp, "report.html")
    _remove_html_element(report_fp, "p", elem_id="icarus")

    # Remove hyperlinks to krona plots in quast HTML report
    _remove_html_element(report_fp, "p", elem_id="krona")

    # remove "Main menu" button from contig browser
    # (otherwise, the visualization can get messed up once clicked)
    contig_browser_fp = os.path.join(
        results_dp, "icarus_viewers", "contig_size_viewer.html"
    )
    if os.path.exists(contig_browser_fp):
        _remove_html_element(contig_browser_fp, "div", elem_id="to_main_menu_button")

    # make all the external links open in a new tab
    _modify_links(report_fp)


def _zip_dir(zip_object: ZipFile, directory: str) -> None:
    for root, _, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            arcname = os.path.relpath(file_path, os.path.dirname(directory))
            zip_object.write(file_path, arcname=arcname)


def _zip_additional_reports(path_to_dirs: list, output_filename: str) -> None:
    with ZipFile(output_filename, "w") as zipf:
        for directory in path_to_dirs:
            _zip_dir(zipf, directory)


def _visualize_quast(
    output_dir: str,
    contigs: ContigSequencesDirFmt,
    min_contig: int = 500,
    threads: int = 1,
    k_mer_stats: bool = False,
    k_mer_size: int = 101,
    contig_thresholds: List[int] = [0, 1000, 5000, 10000, 250000, 500000],
    memory_efficient: bool = False,
    min_alignment: int = 65,
    min_identity: float = 90.0,
    ambiguity_usage: str = "one",
    ambiguity_score: float = 0.99,
    no_icarus: bool = False,
    reads: Union[
        SingleLanePerSamplePairedEndFastqDirFmt, SingleLanePerSampleSingleEndFastqDirFmt
    ] = None,
    references: DNAFASTAFormat = None,
    mapped_reads: BAMDirFmt = None,
) -> None:
    kwargs = {
        k: v
        for k, v in locals().items()
        if k
        not in ["output_dir", "contigs", "reads", "references", "mapped_reads", "ctx"]
    }

    common_args = _process_common_input_params(
        processing_func=_process_quast_arg, params=kwargs
    )

    reads_fps = {}
    paired = False
    if reads:
        paired = isinstance(reads, SingleLanePerSamplePairedEndFastqDirFmt)
        manifest = reads.manifest.view(pd.DataFrame)
        for samp in list(manifest.index):
            reads_fps[samp] = {
                "fwd": manifest.loc[samp, "forward"],
                "rev": manifest.loc[samp, "reverse"] if paired else None,
            }

    with tempfile.TemporaryDirectory() as tmp:
        results_dir = os.path.join(tmp, "results")

        # run quast
        samples = _evaluate_contigs(
            results_dir,
            contigs,
            reads_fps,
            paired,
            references,
            mapped_reads,
            common_args,
        )

        tabular_results = _create_tabular_results(results_dir, contig_thresholds)
        tabular_results.to_csv(os.path.join(results_dir, "quast_results.tsv"), sep="\t")

        # fix/remove some URLs
        _fix_html_reports(results_dir)

        # Copy templates to output dir
        copy_tree(os.path.join(TEMPLATES, "quast"), output_dir)

        # Copy results to output dir
        copy_tree(results_dir, os.path.join(output_dir, "quast_data"))

        # Zip summary, not_aligned and runs_per_reference dirs for download
        dirnames = ["not_aligned", "runs_per_reference", "summary"]
        zip_these_dirs = [
            os.path.join(output_dir, "quast_data", f"{dirname}") for dirname in dirnames
        ]
        output_filename = os.path.join(
            output_dir, "quast_data", "additional_reports.zip"
        )
        _zip_additional_reports(zip_these_dirs, output_filename)

        context = {
            "tabs": [
                {"title": "QC report", "url": "index.html"},
            ],
            "samples": json.dumps(samples),
        }

        templates = [
            os.path.join(TEMPLATES, "quast", "index.html"),
        ]
        if not no_icarus:
            templates.append(os.path.join(TEMPLATES, "quast", "q2_icarus.html"))
            context["tabs"].append({"title": "Contig browser", "url": "q2_icarus.html"})
        if os.path.isdir(os.path.join(output_dir, "quast_data", "krona_charts")):
            templates.append(os.path.join(TEMPLATES, "quast", "q2_krona_charts.html"))
            context["tabs"].append(
                {"title": "Krona charts", "url": "q2_krona_charts.html"}
            )

        q2templates.render(templates, output_dir, context=context)


def _create_tabular_results(results_dir: str, contig_thresholds: list) -> pd.DataFrame:
    """
    This function will create the tabular results after QUAST has run.

    Args:
        - results_dir(str): The directory were the results of QUAST are saved.
        - contig_thresholds(list): list of contig thresholds

    Returns:
        a Pandas dataframe with the tabular data.
    """
    subdir = os.path.join(results_dir, "combined_reference")
    if os.path.isdir(subdir):
        report_fp = os.path.join(subdir, "transposed_report.tsv")
    else:
        report_fp = os.path.join(results_dir, "transposed_report.tsv")

    transposed_report = pd.read_csv(report_fp, sep="\t", header=0)
    transposed_report_parsed = _parse_columns(transposed_report, contig_thresholds)
    return transposed_report_parsed


def evaluate_contigs(
    # TODO: expose more parameters
    ctx,
    contigs,
    reads=None,
    references=None,
    mapped_reads=None,
    min_contig=500,
    threads=1,
    k_mer_stats=False,
    k_mer_size=101,
    contig_thresholds=[0, 1000, 5000, 10000, 25000, 50000],
    memory_efficient=False,
    min_alignment=65,
    min_identity=90.0,
    no_icarus=False,
    ambiguity_usage="one",
    ambiguity_score=0.99,
):
    kwargs = {k: v for k, v in locals().items() if k not in ["contigs", "ctx"]}
    with tempfile.TemporaryDirectory() as tmp:
        # 1. generate the visualization
        _visualize_quast = ctx.get_action("assembly", "_visualize_quast")
        (visualization,) = _visualize_quast(contigs, **kwargs)

        # 2. after the visualization is generated we need to export the files
        # to get the results table out
        visualization_files_path = os.path.join(tmp, "vis_files")
        visualization.export_data(visualization_files_path)
        report_path = os.path.join(
            visualization_files_path, "quast_data", "quast_results.tsv"
        )
        report_df = pd.read_csv(report_path, sep="\t", header=0)

        # 3. read it as a pandas dataframe then we create the QUASTResults
        tabular_results = ctx.make_artifact("QUASTResults", report_df)

    return tabular_results, visualization
