# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
from pathlib import Path
from typing import List, Union

from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt
from skbio import io

from q2_assembly._utils import _construct_param, _get_sample_from_path


def _process_bowtie2build_arg(arg_key, arg_val):
    """Creates a list with argument and its value to be consumed by
        bowtie2-build script.

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
    elif not isinstance(arg_val, list):
        return [_construct_param(arg_key), str(arg_val)]
    else:
        raise NotImplementedError(
            f'Parsing arguments of type "{type(arg_val)}" is not supported.'
        )


def _construct_function_param_value(arg_key: str, arg_val: str):
    """Validates and constructs the 'function' parameters
        consumed by bowtie2.

    Args:
        arg_key (str): Argument name.
        arg_val (str): Argument value. Should be a comma-separated list of
            exactly 3 elements where first element indicates function
            type and the other two stand for function parameters,
            e.g.: "L,1,-0.5".

    Returns:
        validated_function (str): Validated function to be passed to bowtie2.
    """
    param_split = [x.strip() for x in arg_val.split(",")]
    if len(param_split) != 3:
        raise Exception(
            f'Invalid number of elements in function definition of "{arg_key}".'
        )
    elif param_split[0] not in "CLSG":
        raise Exception(
            f'Invalid function type in "{arg_key}": '
            f'{param_split[0]} was given but only "CLSG" '
            f"are allowed."
        )
    return ",".join(param_split)


def _construct_double_list_param_value(arg_key, arg_val):
    """Validates and constructs the two-integer parameters passed to bowtie2.

    Args:
        arg_key (str): Argument name.
        arg_val (str): Argument value. Should be a comma-separated list of
            exactly 2 elements both elements are integers, e.g.: "10,20".

    Returns:
        validated_function (str): Validated string to be passed to bowtie2.
    """
    param_split = [x.strip() for x in arg_val.split(",")]
    if len(param_split) != 2:
        raise Exception(f'Invalid number of elements for "{arg_key}".')
    try:
        return ",".join([str(int(x)) for x in param_split])
    except ValueError:
        raise Exception(
            f'Both values of "{arg_key}" parameter should be '
            f"integers. Provided values were: {arg_val}."
        )


def _process_bowtie2_arg(arg_key, arg_val):
    """Creates a list with argument and its value to be consumed by bowtie2.

    Args:
        arg_key (str): Argument name.
        arg_val: Argument value.

    Returns:
        [converted_arg, arg_value]: List containing a prepared command line
            parameter and its value.
    """
    if arg_key in "drnk" and arg_val:
        arg_key = arg_key if arg_key == "k" else arg_key.capitalize()
        return [f"-{arg_key}", str(arg_val)]
    elif arg_key == "len":
        return ["-L", str(arg_val)]
    # note: `mp` appears to no longer be a list param value, now an int
    # (in version 2.4.4)
    elif arg_key in ["rdg", "rfg"]:
        return [
            _construct_param(arg_key),
            _construct_double_list_param_value(arg_key, arg_val),
        ]
    elif arg_key == "i":
        return ["-i", _construct_function_param_value(arg_key, arg_val)]
    elif arg_key in ["n_ceil", "score_min"]:
        return [
            _construct_param(arg_key),
            _construct_function_param_value(arg_key, arg_val),
        ]
    elif arg_key == "valid_mate_orientations":
        return [_construct_param(arg_val)]
    elif isinstance(arg_val, bool) and arg_val:
        return ["-a"] if arg_key == "a" else [_construct_param(arg_key)]
    elif not isinstance(arg_val, list):
        return [_construct_param(arg_key), str(arg_val)]
    else:
        raise NotImplementedError(
            f'Parsing arguments of type "{type(arg_val)}" is not supported.'
        )


def _get_subdir_from_path(fp: str, input_type: str = "contigs"):
    """Constructs subdir to be created dependent on the input.

    Args:
        fp (str): Path to the original input file.
        input_type (str): Type of input sequences. Can be mags or contigs.

    Returns:
        subdir (str): Subdir to be created, based on the given input.
    """
    if input_type.lower() == "contigs":
        return _get_sample_from_path(fp)
    elif input_type.lower() == "mags":
        return os.path.splitext(fp)[0].split("/")[-2]
    elif input_type.lower() == "mags-derep":
        return ""
    else:
        raise NotImplementedError(f'Input type "{input_type}" is not supported.')


def _assert_inputs_not_empty(fasta_fps: list):
    empty_files = []
    for fp in fasta_fps:
        if not os.path.getsize(fp):
            empty_files.append(Path(fp).name)
    if empty_files:
        msg = (
            f"The following input files were empty: {empty_files}. "
            "Please filter these files from your input and try again."
        )
        raise ValueError(msg)


def _merge_mags_helper(mags: dict, merged_fp: str):
    """
    Merge multiple MAGs from a collection into a single FASTA file.

    This function iterates over a dictionary of MAGs obtained from
    `mags.feature_dict()`, where each MAG is identified by a unique ID
    and associated with a file path. Each MAG's sequences are read,
    modified to prepend the MAG ID to the sequence ID, and then written
    to a single merged FASTA file.

    Args:
        mags: An object that contains a collection of MAGs, which must provide a
            `feature_dict()` method that returns a dictionary mapping MAG IDs to
            their corresponding file paths.
        merged_fp (str): The file path where the merged FASTA file will be written.

    Side Effects:
        Writes to a file specified by `merged_fp`. This file will contain all the
        sequences from the input MAGs, with each sequence ID modified to include
        its source MAG ID.
    """
    with io.open(merged_fp, "w") as merged_f:
        for mag_id, mag_fp in mags.items():
            for seq in io.read(str(mag_fp), format="fasta"):
                seq.metadata["id"] = f"{mag_id}_{seq.metadata['id']}"
                io.write(seq, into=merged_f, format="fasta")


def _merge_mags(
    mags: Union[MAGSequencesDirFmt, MultiMAGSequencesDirFmt], result_dir: str
) -> List[str]:
    """
    Merge multiple MAG sequences into a single FASTA file.

    This function iterates over all the MAGs provided, reads each sequence,
    modifies its ID to include the MAG ID, and writes the sequence to a new
    merged FASTA file. Depending on whether dereplicated MAGs or MAGs per
    sample were provided, a list with a single file path or multi file paths
    will be returned.

    Args:
        mags (Union[MAGSequencesDirFmt, MultiMAGSequencesDirFmt]): The input MAGs.
        result_dir (str): The directory where the merged MAGs will be saved.

    Returns:
        List(str): The list of file path(s) of the merged FASTA file(s).
    """
    if isinstance(mags, MAGSequencesDirFmt):
        merged_fp = os.path.join(result_dir, "merged.fasta")
        _merge_mags_helper(mags.feature_dict(), merged_fp)
        return [merged_fp]
    elif isinstance(mags, MultiMAGSequencesDirFmt):
        all_fps = []
        for sample_id, mags_dict in mags.sample_dict().items():
            sample_dir = os.path.join(result_dir, sample_id)
            os.makedirs(sample_dir, exist_ok=True)
            merged_fp = os.path.join(sample_dir, "merged.fasta")
            _merge_mags_helper(mags_dict, merged_fp)
            all_fps.append(merged_fp)
        return all_fps


def _is_flat_dir(directory: str) -> bool:
    """
    Check if the specified directory contains any subdirectories.

    Args::
        directory (str): The path to the directory to check.

    Returns:
        bool: True if the directory is flat (no subdirectories), False otherwise.
    """
    for entry in os.listdir(directory):
        entry_path = os.path.join(directory, entry)
        if os.path.isdir(entry_path):
            return False
    return True
