# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
import shutil
import subprocess

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


def generate_abundances(
        profile, genomes, sample_name, mu=0, sigma=1, lambd=0.5, random_seed=42
) -> pd.DataFrame:
    np.random.seed(random_seed)
    genome_ids = sorted(list(genomes.file_dict().keys()))
    num_genomes = len(genome_ids)

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

    return pd.DataFrame(
        data={sample_name: abundances}, index=pd.Index(genome_ids, name="id")
    )


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
    sample, genomes, abundances: pd.DataFrame, total_reads, results_dir, threads, read_len, seed
):
    genomes_dict = genomes.file_dict() # {genome1: path1, genome2: path2}
    abundances_dict = abundances[sample].to_dict() # {genome1: abundance1, genome2: abundance2}
    genomes_final = {}
    for genome_id in genomes_dict:
        genomes_final[genome_id] = {
            "fp": genomes_dict[genome_id],
            "reads": int(total_reads * abundances_dict[genome_id]),
        }
    for genome_id, genome_data in genomes_final.items():
        reads_count = genome_data["reads"]
        genome_fp = genome_data["fp"]
        print(f"[Sample: {sample}] Generating {reads_count} reads for {genome_id}")

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
            genome_fp,
            "--num-fragments",
            str(reads_count),
            "--out",
            os.path.join(results_dir, f"{sample}_{genome_id}_00_L001_R1_001.fastq.gz"),
            "--out-right",
            os.path.join(results_dir, f"{sample}_{genome_id}_00_L001_R2_001.fastq.gz"),
        ]

        run_command(cmd, verbose=True)

    # combine all reads into a single file and clean up
    _combine_reads(sample, results_dir, orientation="forward")
    _combine_reads(sample, results_dir, orientation="reverse")


def _simulate_reads_mason(
    reference_genomes: GenomeSequencesDirectoryFormat,
    sample_name: str,
    abundance_profile: str = None,
    abundances: pd.DataFrame = None,
    num_reads: int = 1000000,
    read_length: int = 100,
    random_seed: int = 42,
    threads: int = 1,
) -> (CasavaOneEightSingleLanePerSampleDirFmt, pd.DataFrame):

    # we will copy the genomes to a temporary directory: mason generates fai
    # files that would interfere with validation
    tmp_refs = GenomeSequencesDirectoryFormat()
    shutil.copytree(str(reference_genomes.path), str(tmp_refs.path), dirs_exist_ok=True)

    result_reads = CasavaOneEightSingleLanePerSampleDirFmt()
    abundances = generate_abundances(
        abundance_profile, tmp_refs, sample_name, random_seed=random_seed
    )
    _process_sample(
        sample=sample_name,
        genomes=tmp_refs,
        abundances=abundances,
        total_reads=num_reads,
        results_dir=str(result_reads),
        threads=threads,
        read_len=read_length,
        seed=random_seed,
    )

    del tmp_refs

    return result_reads, abundances.T


def simulate_reads_mason(
    ctx,
    reference_genomes,
    sample_names,
    abundance_profiles=None,
    abundances=None,
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
    merge_tables = ctx.get_action("feature_table", "merge")

    samples, abundances = [], []
    for sample_name, abundance_profile, reads, length in zip(
        sample_names, abundance_profiles, num_reads, read_length
    ):
        reads, abundance = _simulate(
            reference_genomes=reference_genomes,
            abundance_profile=abundance_profile,
            sample_name=sample_name,
            num_reads=reads,
            read_length=length,
            **kwargs,
        )
        samples.append(reads)
        abundances.append(abundance)

    (collated_samples,) = collate_reads(samples)
    (merged_abundances,) = merge_tables(abundances)

    return collated_samples, merged_abundances
