# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import subprocess
import tempfile
import uuid
from typing import List

import pkg_resources
import shortuuid
from bs4 import BeautifulSoup as BS
from skbio import io

EXTERNAL_CMD_WARNING = (
    "Running external command line application(s). "
    "This may print messages to stdout and/or stderr.\n"
    "The command(s) being run are below. These commands "
    "cannot be manually re-run as they will depend on "
    "temporary files that no longer exist."
)


def run_command(cmd, verbose=True, concat=False):
    if verbose:
        print(EXTERNAL_CMD_WARNING)
        print("\nCommand:", end=" ")
        print(" ".join(cmd), end="\n\n")
    if concat:
        subprocess.run(cmd, shell=True, check=True)
    else:
        subprocess.run(cmd, check=True)


def run_commands_with_pipe(cmd1, cmd2, verbose=True):
    """Runs two consecutive commands using a pipe"""
    if verbose:
        print(EXTERNAL_CMD_WARNING)
        print("\nCommand:", end=" ")
        print(f'{" ".join(cmd1)} | {" ".join(cmd2)}', end="\n\n")
    out1 = subprocess.run(cmd1, check=True, capture_output=True)
    subprocess.run(cmd2, input=out1.stdout, check=True)


def _construct_param(arg_name):
    """Converts argument name into a command line parameter."""
    return f'--{arg_name.replace("_", "-")}'


def _process_common_input_params(processing_func, params: dict) -> List[str]:
    """Converts provided arguments and their values.

    Conversion is entirely dependent on the passed 'processing_func'
    that processes individual arguments. The output is a list of
    parameters with their values that can be directly passed to the
    respective command.

    Arguments without any value are skipped.
    Arguments of boolean type are only appended if True.
    Any other argument is processed using the 'processing_func' and
    appended to the final list.

    Args:
        processing_func: Function to be used for formatting a single argument.
        params (dict): Dictionary of parameter: value pairs to be processed.

    Returns:
        processed_args (list): List of processed arguments and their values.

    """
    processed_args = []
    for arg_key, arg_val in params.items():
        if not arg_val:
            continue
        else:
            processed_args.extend(processing_func(arg_key, arg_val))
    return processed_args


def _remove_html_element(fp: str, tag: str, elem_id: str):
    """Removes an HTML tag and its contents.

    Uses BeautifulSoup to open an HTML file, find an element by tag and
    its id and remove it, if found. The original file will be overwritten
    by its modified version.

    Args:
         fp (str): Path to the original HTML file.
         tag (str): Type of the tag to be removed.
         elem_id (str): ID of the element (tag) to be removed.
    """
    with open(fp, "r") as r:
        soup = BS(r.read(), "html.parser")
        element = soup.find(tag, id=elem_id)
        if element:
            element.decompose()
    os.remove(fp)

    with open(fp, "w") as r:
        r.write(str(soup))


def _modify_links(fp: str):
    """Modifies all "a" tags to automatically open in a new tab rather
        than the original iFrame.

    Uses BeautifulSoup to open an HTML file, find all "a" tags and
    add a "target" property. The original file will be overwritten
    by its modified version.

    Args:
         fp (str): Path to the original HTML file.
    """
    with open(fp, "r") as r:
        soup = BS(r.read(), "html.parser")
        links = soup.find_all("a")
        for line in links:
            line["target"] = "_blank"
    os.remove(fp)

    with open(fp, "w") as r:
        r.write(str(soup))


def _get_sample_from_path(fp):
    """Extracts sample name from a contig's file path."""
    return os.path.basename(fp).rsplit("_contigs.fa", maxsplit=1)[0]


def get_relative_data_path(package, filename):
    """Get data path relative to the provided package.

    Args:
        package (str): The package we are getting the data path under
        filename (str): The name of the file/dir we are trying to get
    """
    return pkg_resources.resource_filename(package, f"data/{filename}")


def get_file_extension(filepath):
    """Extract the extension of a file to see if it is compressed
    or not.

    Args:
        filepath (str): the path of the file to get the extension

    """
    ext = ""
    parts = filepath.split(".")

    if parts[-1] == "gz":
        ext += f".{parts[-2]}.gz"
    else:
        # take only the last part that
        # will correspond to the extension
        ext = f".{parts[-1]}"
    return ext


def concatenate_files(input_files, output_file):
    """Concatenate the content of the files in input_files and
    save the content in the output_file.

    Args:
        input_files (list): list of all files to be concatenated
        output_file (str): the path to the resulting file
    """
    cmd = ["cat", *input_files]
    subprocess.run(cmd, stdout=open(output_file, "w"), check=True)


def generate_shortuuid():
    """
    Generate a shortuuid string.

    Returns: a shortuuid strings.
    """
    return shortuuid.uuid()


def generate_uuid3(namespace, name):
    """
    Generate a UUID3 string.
    Args:
        - namespace: UUID namespace.
        - name: Name to be used in the UUID generation.
    Returns: UUID3 string.
    """
    return uuid.uuid3(namespace, name)


def generate_uuid4():
    """
    Generate a UUID4 string.

    Returns: UUID4 string.
    """
    return uuid.uuid4()


def generate_uuid5(namespace, name):
    """
    Generate a UUID5 string.
    Args:
        - namespace: UUID namespace.
        - name: Name to be used in the UUID generation.
    Returns: UUID5 string.
    """
    return uuid.uuid5(namespace, name)


def _generate_unique_uuid(uuid_func, namespace, contig_id, new_ids, uuid_type):
    """Generates a unique UUID for a contig.
    Args:
        - uuid_func: UUID generation function.
        - namespace: UUID namespace.
        - contig_id: Contig ID, what will be changed.
        - new_ids: Set of already generated UUIDs.
        - uuid_type: the type of uuid choose by the user

    Returns: the new generated contig_id
    """
    while True:
        if uuid_type in ["uuid4", "shortuuid"]:
            new_id = str(uuid_func())
        else:
            new_id = str(uuid_func(namespace, contig_id))

        if new_id not in new_ids:
            new_ids.add(new_id)
            return new_id


def modify_contig_ids(out, sample, uuid_type: str):
    """Modifies the contig IDs to include the sample name and UUID.

    Args:
        out: output format
        sample: Name of the sample to be processed.
        uuid_type: UUID type to be used in the contig ID.
    """
    uuid_funcs = {
        "uuid3": generate_uuid3,
        "uuid4": generate_uuid4,
        "uuid5": generate_uuid5,
        "shortuuid": generate_shortuuid,
    }
    uuid_func = uuid_funcs.get(uuid_type)

    path_to_contigs = os.path.join(str(out), f"{sample}_contigs.fa")
    with tempfile.TemporaryDirectory() as tmp:
        path_to_contigs_mod = os.path.join(tmp, f"{sample}_modified_contigs.fa")
        namespace = uuid.NAMESPACE_OID
        new_ids = set()
        with io.open(path_to_contigs_mod, "w") as modified_contigs:
            for contig in io.read(path_to_contigs, format="fasta"):
                contig.metadata["id"] = _generate_unique_uuid(
                    uuid_func, namespace, contig.metadata["id"], new_ids, uuid_type
                )
                io.write(contig, format="fasta", into=modified_contigs)

        all_contigs = io.read(path_to_contigs_mod, format="fasta")
        with io.open(path_to_contigs, "w") as contigs_file:
            io.write(all_contigs, format="fasta", into=contigs_file)
