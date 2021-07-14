# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os

from q2_assembly._utils import _construct_param


def _process_bowtie2_arg(arg_key, arg_val):
    """Creates a list with argument and its value.

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
            f'Parsing arguments of type "{type(arg_val)}" is not supported.')


def _get_subdir_from_path(fp: str, input_type: str = 'contigs'):
    """Constructs subdir to be created dependent on the input.

    Args:
        fp (str): Path to the original input file.
        input_type (str): Type of input sequences. Can be mags or contigs.

    Returns:
        subdir (str): Subdir to be created, based on the given input.
    """
    if input_type.lower() == 'contigs':
        return os.path.basename(fp).split('_', maxsplit=1)[0]
    elif input_type.lower() == 'mags':
        fpl = os.path.splitext(fp)
        return os.path.join(*fpl[0].split('/')[-2:])
    else:
        raise NotImplementedError(f'Input type "{input_type}" '
                                  f'is not supported.')
