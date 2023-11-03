# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import os
import shutil
import subprocess
import tempfile
from copy import deepcopy
from typing import List

import biom
import pandas as pd
import skbio
from q2_types.feature_data import DNAFASTAFormat, DNAIterator
from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt

from .._utils import _process_common_input_params, run_command


def _process_iss_arg(arg_key, arg_val):
    """Creates a list with argument and its value.

    Argument values represented by a list will be converted to a single
    string joined by spaces, e.g.: [1, 2, 3] -> '1 2 3'.
    Argument names will be converted to command line parameters by
    appending a '--' prefix, e.g.: 'some_parameter' -> '--some_parameter'.

    Args:
        arg_key (str): Argument name.
        arg_val: Argument value.

    Returns:
        [converted_arg, arg_value, ...]: List containing a prepared command
            line parameter and its value(s).
    """
    if isinstance(arg_val, bool) and arg_val:
        return [f"--{arg_key}"]
    elif not isinstance(arg_val, list):
        return [f"--{arg_key}", str(arg_val)]
    else:
        flags = [f"--{arg_key}"]
        flags.extend([str(x) for x in arg_val])
        return flags


def _rename_reads_files(dp):
    reads = sorted(glob.glob(os.path.join(dp, "*.fastq.gz")))
    reads_new = [r.replace(".fastq.gz", "_001.fastq.gz") for r in reads]
    for r, rn in zip(reads, reads_new):
        os.rename(r, rn)
    return reads_new


def _generate_reads(samples, args, result_fp):
    base_cmd = ["iss", "generate", "--compress"]
    base_cmd.extend(args)

    for s in samples:
        cmd = deepcopy(base_cmd)
        sample_prefix = os.path.join(result_fp, f"{s}_00_L001")

        cmd.extend(["--output", sample_prefix])

        try:
            run_command(cmd, verbose=True)
        except subprocess.CalledProcessError as e:
            raise Exception(
                "An error was encountered while running InSilicoSeq, "
                f"(return code {e.returncode}), please inspect "
                "stdout and stderr to learn more."
            )

    _rename_reads_files(result_fp)


def _abundances_to_biom(abundance_fps):
    abundances = []
    for f in abundance_fps:
        sample = os.path.splitext(os.path.basename(f))[0].split("_")[0]
        abundances.append(
            pd.read_csv(f, sep="\t", index_col=0, header=None, names=[sample])
        )
    biom_df = pd.concat(abundances, axis=1)
    return biom.Table(
        data=biom_df.values,
        observation_ids=biom_df.index.tolist(),
        sample_ids=biom_df.columns.tolist(),
    )


def _ensure_sample_names_exists(sample_names):
    if not sample_names:
        # If it's empty or None, create a list with a default element "sample"
        sample_names = ["sample"]
        print(
            'The "--p-sample-names" option was not provided. '
            'Only one sample will be created with the prefix "sample".'
            "\n"
        )
    return sample_names


# TODO: allow to input custom genomes and sample from those
def generate_reads(
    genomes: DNAFASTAFormat = None,
    sample_names: List[str] = None,
    n_genomes: int = 10,
    ncbi: List[str] = ["bacteria"],
    n_genomes_ncbi: List[int] = [10],
    abundance: str = "lognormal",
    coverage: str = "off",
    n_reads: int = 1000000,
    mode: str = "kde",
    model: str = "HiSeq",
    gc_bias: bool = False,
    cpus: int = 1,
    debug: bool = False,
    seed: int = 0,
) -> (CasavaOneEightSingleLanePerSampleDirFmt, DNAFASTAFormat, biom.Table):
    if coverage == "off":
        coverage = None
    _locals = locals().copy()
    available_genomes = 0
    if genomes:
        if n_genomes_ncbi or ncbi:
            print(
                'Template genome sequences were provided - "n-genomes-ncbi" '
                'and "ncbi" parameters will be ignored.'
            )
            _locals["n_genomes_ncbi"], _locals["ncbi"] = None, None
        for _ in genomes.view(DNAIterator):
            available_genomes += 1
    elif n_genomes_ncbi and ncbi and (len(n_genomes_ncbi) != len(ncbi)):
        raise Exception(
            "Genome counts (--n_genomes_ncbi) need to correspond "
            "to the kingdoms names (--ncbi). You provided "
            f"{len(ncbi)} kingdom(s) but {len(n_genomes_ncbi)} "
            f"corresponding genome counts were found. Please "
            f"correct your input."
        )

    if genomes and (n_genomes >= available_genomes):
        print(
            f"The number of available genomes ({available_genomes}) is "
            f"smaller than the requested number of genomes per sample "
            f"({n_genomes}). The number of requested genomes will be "
            f"reduced to {available_genomes - 1}."
        )
        _locals["n_genomes"] = available_genomes - 1

    sample_names = _ensure_sample_names_exists(sample_names)

    if len(set(sample_names)) < len(sample_names):
        dupl = {str(x) for x in sample_names if sample_names.count(x) > 1}
        raise Exception(
            "Sample names need to be unique. Found duplicated "
            f'names: {", ".join(sorted(dupl))}'
        )

    kwargs = {k: v for k, v in _locals.items() if k not in ["sample_names"]}
    args = _process_common_input_params(processing_func=_process_iss_arg, params=kwargs)

    with tempfile.TemporaryDirectory() as tmp:
        result_reads = CasavaOneEightSingleLanePerSampleDirFmt()
        result_genomes = DNAFASTAFormat()

        # simulate reads
        _generate_reads(sample_names, args, tmp)

        # move reads into CasavaFmt
        for f in glob.glob(os.path.join(tmp, "*.fastq.gz")):
            shutil.move(f, str(result_reads))

        # move original genomes into DNAFASTAFmt if found
        # otherwise return empty file
        genome_ids = []
        with result_genomes.open() as fout:
            for f in sorted(glob.glob(os.path.join(tmp, "*.fasta"))):
                for seq in skbio.read(f, format="fasta"):
                    if seq.metadata["id"] not in genome_ids:
                        genome_ids.append(seq.metadata["id"])
                        seq.write(fout)

        # convert abundances to a biom table
        abund_suffix = "coverage" if coverage else "abundance"
        abundance_fps = sorted(glob.glob(os.path.join(tmp, f"*_{abund_suffix}.txt")))
        result_biom = _abundances_to_biom(abundance_fps)

        return result_reads, result_genomes, result_biom
