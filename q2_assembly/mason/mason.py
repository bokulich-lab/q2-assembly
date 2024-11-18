# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import shutil
import subprocess
import tempfile
from typing import List

from q2_types.genome_data import GenomeData
from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt

from .._utils import run_command


def _process_mason_arg(arg_key, arg_val):
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


def _simulate_reads(samples, args, result_fp):
    base_cmd = ["mason_simulator"]
    base_cmd.extend(args)

    for s in samples:
        cmd = base_cmd.copy()
        sample_prefix = os.path.join(result_fp, f"{s}_00_L001")

        cmd.extend(["--output", sample_prefix])

        try:
            run_command(cmd, verbose=True)
        except subprocess.CalledProcessError as e:
            raise Exception(
                "An error was encountered while running Mason, "
                f"(return code {e.returncode}), please inspect "
                "stdout and stderr to learn more."
            )


def simulate_reads(
    reference_genomes: GenomeData,
    sample_names: List[str],
    num_reads: int = 1000000,
    read_length: int = 100,
    fragment_mean_size: int = 500,
    fragment_size_stddev: int = 50,
    error_rate: float = 0.01,
    random_seed: int = 42,
    haplotype_count: int = 1,
    haplotype_diversity: float = 0.0,
    indel_probability: float = 0.0,
    indel_max_length: int = 0,
    snp_probability: float = 0.0,
    snp_transition_transversion_ratio: float = 2.0,
) -> CasavaOneEightSingleLanePerSampleDirFmt:
    sample_names = sample_names or ["sample"]

    if len(set(sample_names)) < len(sample_names):
        dupl = {str(x) for x in sample_names if sample_names.count(x) > 1}
        raise Exception(
            "Sample names need to be unique. Found duplicated "
            f'names: {", ".join(sorted(dupl))}'
        )

    args = _process_mason_arg("num_reads", num_reads)
    args.extend(_process_mason_arg("read_length", read_length))
    args.extend(_process_mason_arg("fragment_mean_size", fragment_mean_size))
    args.extend(_process_mason_arg("fragment_size_stddev", fragment_size_stddev))
    args.extend(_process_mason_arg("error_rate", error_rate))
    args.extend(_process_mason_arg("random_seed", random_seed))
    args.extend(_process_mason_arg("haplotype_count", haplotype_count))
    args.extend(_process_mason_arg("haplotype_diversity", haplotype_diversity))
    args.extend(_process_mason_arg("indel_probability", indel_probability))
    args.extend(_process_mason_arg("indel_max_length", indel_max_length))
    args.extend(_process_mason_arg("snp_probability", snp_probability))
    args.extend(_process_mason_arg("snp_transition_transversion_ratio", snp_transition_transversion_ratio))

    with tempfile.TemporaryDirectory() as tmp:
        result_reads = CasavaOneEightSingleLanePerSampleDirFmt()

        # simulate reads
        for genome in reference_genomes.file.view(DNAIterator):
            genome_fp = os.path.join(tmp, f"{genome.metadata['id']}.fasta")
            with open(genome_fp, 'w') as f:
                genome.write(f)
            args.extend(["--input-reference", genome_fp])
            _simulate_reads(sample_names, args, tmp)

        # move reads into CasavaFmt
        for f in os.listdir(tmp):
            shutil.move(os.path.join(tmp, f), str(result_reads))

        return result_reads
