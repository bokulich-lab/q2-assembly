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
from typing import List

import numpy as np
import pandas as pd
from q2_types.genome_data import GenomeSequencesDirectoryFormat
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


def generate_abundances(profile, num_genomes, mu=0, sigma=1, lambd=0.5, random_seed=42):
    np.random.seed(random_seed)

    if profile == "uniform":
        abundances = [1 / num_genomes] * num_genomes
    elif profile == "exponential":
        abundances = np.exp(-lambd * np.arange(num_genomes))
        abundances /= abundances.sum()
    elif profile == "lognormal":
        abundances = np.exp(mu + sigma * np.random.randn(num_genomes))
        abundances /= abundances.sum()
    else:
        raise ValueError(f"Invalid abundance profile option: '{profile}'")

    return list(abundances)


def abundances_to_df(abundances, genome_files, sample_id):
    data = {
        "id": [os.path.basename(f).split(".")[0] for f in genome_files],
        sample_id: abundances,
    }
    df = pd.DataFrame(data)
    df.set_index("id", inplace=True)
    return df


def _combine_reads(sample_id, results_dir, orientation="forward"):
    _idx = "1" if orientation == "forward" else "2"

    all_files = sorted(
        glob.glob(
            os.path.join(results_dir, f"{sample_id}_*_00_L001_R{_idx}_001.fastq.gz")
        )
    )
    with open(
        os.path.join(results_dir, f"{sample_id}_00_L001_R{_idx}_001.fastq.gz"), "wb"
    ) as out:
        cmd = ["cat", *all_files]
        subprocess.run(cmd, check=True, stdout=out)

    # clean up
    for f in all_files:
        os.remove(f)


def _process_sample(
    sample, genome_files, abundances, total_reads, results_dir, threads, read_len, seed
):
    for genome_file, abundance in zip(genome_files, abundances):
        genome_reads = int(total_reads * abundance)
        _id = os.path.splitext(os.path.basename(genome_file))[0]
        print(f"[Sample: {sample}] Generating {genome_reads} reads for {_id}")

        cmd = [
            "mason_simulator",
            "-v",
            "--seed",
            str(seed),
            "--illumina-read-length",
            str(read_len),
            "--seq-technology",
            "illumina",
            "--num-threads",
            str(threads),
            "--input-reference",
            genome_file,
            "--num-fragments",
            str(genome_reads),
            "--out",
            os.path.join(results_dir, f"{sample}_{_id}_00_L001_R1_001.fastq.gz"),
            "--out-right",
            os.path.join(results_dir, f"{sample}_{_id}_00_L001_R2_001.fastq.gz"),
        ]

        run_command(cmd, verbose=True)

    # combine all reads into a single file and clean up
    _combine_reads(sample, results_dir, orientation="forward")
    _combine_reads(sample, results_dir, orientation="reverse")


def _simulate_reads_mason(
    reference_genomes: GenomeSequencesDirectoryFormat,
    sample_name: str,
    abundance_profile: str = None,
    num_reads: int = 1000000,
    read_length: int = 100,
    random_seed: int = 42,
    threads: int = 1,
) -> CasavaOneEightSingleLanePerSampleDirFmt:

    # we will copy the genomes to a temporary directory: mason generates fai
    # files that would interfere with validation
    tmp_refs = GenomeSequencesDirectoryFormat()
    shutil.copytree(str(reference_genomes.path), str(tmp_refs.path), dirs_exist_ok=True)

    result_reads = CasavaOneEightSingleLanePerSampleDirFmt()
    genome_files = list(tmp_refs.file_dict().values())
    abundances = generate_abundances(
        abundance_profile, len(genome_files), random_seed=random_seed
    )
    _process_sample(
        sample=sample_name,
        genome_files=genome_files,
        abundances=abundances,
        total_reads=num_reads,
        results_dir=str(result_reads),
        threads=threads,
        read_len=read_length,
        seed=random_seed,
    )

    del tmp_refs

    return result_reads


def simulate_reads_mason(
    ctx,
    reference_genomes,
    sample_names,
    abundance_profiles,
    num_reads=[1000000],
    read_length=[100],
    random_seed=42,
    threads=1,
    num_partitions=None,
):
    kwargs = {
        k: v
        for k, v in locals().items()
        if k
        not in [
            "reference_genomes",
            "sample_names",
            "num_partitions",
            "ctx",
            "abundance_profiles",
            "num_reads",
            "read_length",
        ]
    }

    sample_names = sample_names or ["sample"]
    if len(set(sample_names)) < len(sample_names):
        dupl = {str(x) for x in sample_names if sample_names.count(x) > 1}
        raise ValueError(
            "Sample names need to be unique. Found duplicated "
            f'names: {", ".join(sorted(dupl))}'
        )

    if len(abundance_profiles) != len(sample_names) and len(abundance_profiles) != 1:
        raise ValueError(
            f"The length of abundance_profiles must be either 1 or equal to the "
            f"number of sample names. Provided: {len(abundance_profiles)}, "
            f"Expected: 1 or {len(sample_names)}"
        )

    if len(num_reads) != len(sample_names) and len(num_reads) != 1:
        raise ValueError(
            f"The length of num_reads must be either 1 or equal to the "
            f"number of sample names. Provided: {len(num_reads)}, "
            f"Expected: 1 or {len(sample_names)}"
        )

    if len(read_length) != len(sample_names) and len(read_length) != 1:
        raise ValueError(
            f"The length of read_length must be either 1 or equal to the "
            f"number of sample names. Provided: {len(read_length)}, "
            f"Expected: 1 or {len(sample_names)}"
        )

    if len(abundance_profiles) == 1:
        abundance_profiles = abundance_profiles * len(sample_names)

    if len(num_reads) == 1:
        num_reads = num_reads * len(sample_names)

    if len(read_length) == 1:
        read_length = read_length * len(sample_names)

    _simulate = ctx.get_action("assembly", "_simulate_reads_mason")
    collate_reads = ctx.get_action("fondue", "combine_seqs")

    samples = []
    for sample_name, abundance_profile, reads, length in zip(
        sample_names, abundance_profiles, num_reads, read_length
    ):
        (sample,) = _simulate(
            reference_genomes=reference_genomes,
            abundance_profile=abundance_profile,
            sample_name=sample_name,
            num_reads=reads,
            read_length=length,
            **kwargs,
        )
        samples.append(sample)

    (collated_samples,) = collate_reads(samples)

    return collated_samples
