# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os

from q2_assembly._utils import _construct_param


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
            "Invalid number of elements in function definition " f'of "{arg_key}".'
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
    elif arg_key in ["mp", "rdg", "rfg"]:
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
        return os.path.basename(fp.replace('_contigs.fa', ''))
    elif input_type.lower() == "mags":
        fpl = os.path.splitext(fp)
        return os.path.join(*fpl[0].split("/")[-2:])
    else:
        raise NotImplementedError(f'Input type "{input_type}" ' f"is not supported.")
