# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import os
import subprocess
from copy import deepcopy

from q2_types.bowtie2 import Bowtie2IndexDirFmt
from q2_types_genomics.per_sample_data import (ContigSequencesDirFmt,
                                               MultiBowtie2IndexDirFmt,
                                               MultiMAGSequencesDirFmt)

from .._utils import (run_command, _construct_param,
                      _process_common_input_params)


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


def _index_seqs(
        fasta_fps: list, result_fp: str, common_args: list,
        input_type: str = 'contigs'
):
    """Runs the indexing using bowtie2

    Constructs and runs the final bowtie2-build command.

    Args:
        fasta_fps (list): List of FASTA files to be indexed.
        result_fp (str): Path to the result file where the indices
            will be created.
        common_args (list): List of common flags and their values for
            the bowtie2-build command.
        input_type (str): Type of input sequences. Can be mags or contigs.
    """
    base_cmd = ['bowtie2-build']
    base_cmd.extend(common_args)

    for fp in fasta_fps:
        sample_dp = os.path.join(
            result_fp, _get_subdir_from_path(fp, input_type))
        os.makedirs(sample_dp)

        cmd = deepcopy(base_cmd)
        cmd.extend([fp, os.path.join(sample_dp, 'index')])

        try:
            run_command(cmd, verbose=True)
        except subprocess.CalledProcessError as e:
            raise Exception(
                'An error was encountered while running Bowtie2, '
                f'(return code {e.returncode}), please inspect '
                'stdout and stderr to learn more.')


def index_contigs(
        contigs: ContigSequencesDirFmt,
        large_index: bool = False,
        debug: bool = False,
        sanitized: bool = False,
        verbose: bool = False,
        noauto: bool = False,
        packed: bool = False,
        bmax: int = None,
        bmaxdivn: int = None,
        dcv: int = None,
        nodc: bool = False,
        offrate: int = None,
        ftabchars: int = None,
        threads: int = None,
        seed: int = None
) -> Bowtie2IndexDirFmt:
    kwargs = {k: v for k, v in locals().items() if k not in ['contigs']}
    common_args = _process_common_input_params(
        processing_func=_process_bowtie2_arg, params=kwargs
    )
    result = Bowtie2IndexDirFmt()

    contig_fps = sorted(glob.glob(os.path.join(str(contigs), '*_contigs.fa')))
    _index_seqs(contig_fps, str(result), common_args, 'contigs')

    return result


def index_mags(
        mags: MultiMAGSequencesDirFmt,
        large_index: bool = False,
        debug: bool = False,
        sanitized: bool = False,
        verbose: bool = False,
        noauto: bool = False,
        packed: bool = False,
        bmax: int = None,
        bmaxdivn: int = None,
        dcv: int = None,
        nodc: bool = False,
        offrate: int = None,
        ftabchars: int = None,
        threads: int = None,
        seed: int = None
) -> MultiBowtie2IndexDirFmt:
    kwargs = {k: v for k, v in locals().items() if k not in ['mags']}
    common_args = _process_common_input_params(
        processing_func=_process_bowtie2_arg, params=kwargs
    )

    result = MultiBowtie2IndexDirFmt()

    mag_fps = sorted(glob.glob(os.path.join(str(mags), '*', '*.fa*')))
    _index_seqs(mag_fps, str(result), common_args, 'mags')

    return result
