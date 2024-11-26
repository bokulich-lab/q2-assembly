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


def generate_abundances(profiles, num_genomes, mu=0, sigma=1, lambd=0.5):
    all_abundances = []
    for profile in profiles:
        abundances = []
        if profile == "uniform":
            abundances = [1 / num_genomes] * num_genomes
        elif profile == "exponential":
            abundances = np.exp(-lambd * np.arange(num_genomes))
            abundances /= abundances.sum()
        elif profile == "lognormal":
            abundances = np.exp(mu + sigma * np.random.randn(num_genomes))
            abundances /= abundances.sum()
        else:
            print(f"Invalid abundance profile option: {profile}")
        all_abundances.append(abundances)
    return all_abundances


def abundances_to_df(abundances, genome_files, sample_id):
    data = {
        "id": [os.path.basename(f).split(".")[0] for f in genome_files],
        sample_id: abundances,
    }
    df = pd.DataFrame(data)
    df.set_index("id", inplace=True)
    return df


def _process_sample(
    sample, genome_files, abundances, total_reads, tmp_dir, threads, read_len, seed
):
    for genome_file, abundance in zip(genome_files, abundances):
        genome_reads = int(total_reads * abundance)
        _id = os.path.basename(genome_file).replace(".fasta", "")
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
            os.path.join(tmp_dir, f"{sample}_{_id}_L001_R1_001.fastq.gz"),
            "--out-right",
            os.path.join(tmp_dir, f"{sample}_{_id}_L001_R2_001.fastq.gz"),
        ]

        run_command(cmd, verbose=True)

    # combine all reads into a single file
    cmd = [
        "cat",
        os.path.join(tmp_dir, f"{sample}_*_L001_R1_001.fastq.gz"),
        ">",
        os.path.join(tmp_dir, f"{sample}_L001_R1_001.fastq.gz"),
    ]
    run_command(cmd, verbose=True, concat=True)

    cmd = [
        "cat",
        os.path.join(tmp_dir, f"{sample}_*_L001_R2_001.fastq.gz"),
        ">",
        os.path.join(tmp_dir, f"{sample}_L001_R2_001.fastq.gz"),
    ]
    run_command(cmd, verbose=True, concat=True)


def _simulate_reads_mason(
    reference_genomes: GenomeSequencesDirectoryFormat,
    sample_names: List[str],
    num_reads: int = 1000000,
    read_length: int = 100,
    fragment_mean_size: int = 500,
    fragment_size_stddev: int = 50,
    error_rate: float = 0.01,
    random_seed: int = 42,
    abundance_profiles: List[str] = None,
    threads: int = 1,
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

    with tempfile.TemporaryDirectory() as tmp:
        result_reads = CasavaOneEightSingleLanePerSampleDirFmt()

        genome_files = reference_genomes.genome_dict().values()
        abundances = generate_abundances(abundance_profiles, len(genome_files))
        for sample_name, abundance in zip(sample_names, abundances):
            _process_sample(
                sample_name,
                genome_files,
                abundance,
                num_reads,
                tmp,
                threads,
                read_length,
                random_seed,
            )

        # move reads into CasavaFmt
        for f in os.listdir(tmp):
            shutil.move(os.path.join(tmp, f), str(result_reads))

    for f in glob.glob(os.path.join(str(reference_genomes.path), "*.fai")):
        os.remove(f)

    return result_reads


def simulate_reads_mason(
    ctx,
    reference_genomes,
    sample_names,
    abundance_profiles,
    num_reads=1000000,
    read_length=100,
    fragment_mean_size=500,
    fragment_size_stddev=50,
    error_rate=0.01,
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
        ]
    }
    # check that the length of sample_names and abundance_profiles match
    if len(sample_names) != len(abundance_profiles):
        raise ValueError(
            "The number of sample names and abundance profiles must match."
        )
    _simulate = ctx.get_action("assembly", "_simulate_reads_mason")
    collate_reads = ctx.get_action("fondue", "combine_seqs")

    samples = []
    for sample_name, abundance_profile in zip(sample_names, abundance_profiles):
        kwargs["abundance_profiles"] = [abundance_profile]
        (sample,) = _simulate(
            reference_genomes=reference_genomes,
            sample_names=[sample_name],
            **kwargs,
        )
        samples.append(sample)

    (collated_samples,) = collate_reads(samples)

    return collated_samples
