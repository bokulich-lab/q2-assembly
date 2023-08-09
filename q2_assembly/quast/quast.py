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

import pandas as pd
import pkg_resources
import q2templates
from q2_types.feature_data import DNAFASTAFormat, DNAIterator
from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt

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

    for fp in sorted(glob.glob(os.path.join(str(contigs), "*_contigs.fa"))):
        cmd.append(fp)
        samples.append(_get_sample_from_path(fp))

    if reads:
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

    A few HTML elements generated by QUAST could break the visualisation
    when clicked so they need to be removed or disabled.
    Modifies the files in-place.

    Args:
         results_dp (str): Directory where the reports are stored.

    """

    # remove link to icarus browser from report.html
    # (otherwise, the visualisation can get messed up once clicked)
    report_fp = os.path.join(results_dp, "report.html")
    _remove_html_element(report_fp, "p", elem_id="icarus")

    # remove "Main menu" button from contig browser
    # (otherwise, the visualisation can get messed up once clicked)
    contig_browser_fp = os.path.join(
        results_dp, "icarus_viewers", "contig_size_viewer.html"
    )
    _remove_html_element(contig_browser_fp, "div", elem_id="to_main_menu_button")

    # make all the external links open in a new tab
    _modify_links(report_fp)


def evaluate_contigs(
    # TODO: expose more parameters
    output_dir: str,
    contigs: ContigSequencesDirFmt,
    reads: Union[
        SingleLanePerSamplePairedEndFastqDirFmt, SingleLanePerSampleSingleEndFastqDirFmt
    ] = None,
    references: DNAFASTAFormat = None,
    min_contig: int = 500,
    threads: int = 1,
    k_mer_stats: bool = False,
    k_mer_size: int = 101,
    contig_thresholds: List[int] = [0, 1000, 5000, 10000, 25000, 50000],
):
    kwargs = {
        k: v
        for k, v in locals().items()
        if k not in ["output_dir", "contigs", "reads", "references"]
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
            results_dir, contigs, reads_fps, paired, references, common_args
        )

        # fix/remove some URLs
        _fix_html_reports(results_dir)

        copy_tree(os.path.join(TEMPLATES, "quast"), output_dir)
        copy_tree(results_dir, os.path.join(output_dir, "quast_data"))

        context = {
            "tabs": [
                {"title": "QC report", "url": "index.html"},
                {"title": "Contig browser", "url": "q2_icarus.html"},
            ],
            "samples": json.dumps(samples),
        }

        index = os.path.join(TEMPLATES, "quast", "index.html")
        icarus = os.path.join(TEMPLATES, "quast", "q2_icarus.html")

        templates = [index, icarus]
        q2templates.render(templates, output_dir, context=context)
