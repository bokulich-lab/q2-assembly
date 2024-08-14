# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import shutil
import warnings
from typing import Union

import numpy as np
import skbio.io
from q2_types._util import DNAFASTAFormat
from q2_types.bowtie2 import Bowtie2IndexDirFmt
from q2_types.feature_data import DNAIterator
from q2_types.genome_data import GenomeSequencesDirectoryFormat
from q2_types.per_sample_sequences import BAMDirFmt, ContigSequencesDirFmt
from qiime2.util import duplicate

from q2_assembly._utils import modify_contig_ids


def partition_contigs(
    contigs: ContigSequencesDirFmt, num_partitions: int = None
) -> ContigSequencesDirFmt:
    partitioned_contigs = {}
    contigs = [
        (sample_id, sample_fp) for sample_id, sample_fp in contigs.sample_dict().items()
    ]

    num_samples = len(contigs)
    if num_partitions is None:
        num_partitions = num_samples
    elif num_partitions > num_samples:
        warnings.warn(
            "You have requested a number of partitions"
            f" '{num_partitions}' that is greater than your number"
            f" of samples '{num_samples}.' Your data will be"
            f" partitioned by sample into '{num_samples}'"
            " partitions."
        )
        num_partitions = num_samples

    contigs = np.array_split(contigs, num_partitions)
    for i, samples in enumerate(contigs, 1):
        result = ContigSequencesDirFmt()

        for sample_id, sample_fp in samples:
            duplicate(sample_fp, result.path / os.path.basename(sample_fp))

        # If num_partitions == num_samples we will only have gone through one
        # sample in the above loop and will use its id as a key. Otherwise we
        # may have gone through multiple samples in the above loop and will be
        # using indices for keys
        if num_partitions == num_samples:
            partitioned_contigs[sample_id] = result
        else:
            partitioned_contigs[i] = result

    return partitioned_contigs


def rename_contigs(
    contigs: ContigSequencesDirFmt, uuid_type: str
) -> ContigSequencesDirFmt:
    renamed_contigs = ContigSequencesDirFmt()

    for contigs_fp in contigs.sample_dict().values():
        shutil.copyfile(
            contigs_fp,
            os.path.join(renamed_contigs.path, os.path.basename(contigs_fp)),
        )

    for sample_id, contigs_fp in renamed_contigs.sample_dict().items():
        modify_contig_ids(contigs_fp, sample_id, uuid_type)

    return renamed_contigs


def collate_contigs(contigs: ContigSequencesDirFmt) -> ContigSequencesDirFmt:
    collated_contigs = ContigSequencesDirFmt()

    for contig in contigs:
        for fp in contig.path.iterdir():
            duplicate(fp, collated_contigs.path / fp.name)

    return collated_contigs


def collate_indices(indices: Bowtie2IndexDirFmt) -> Bowtie2IndexDirFmt:
    collated_indices = Bowtie2IndexDirFmt()

    for index in indices:
        for in_sample_dir in index.path.iterdir():
            out_sample_dir = collated_indices.path / in_sample_dir.name
            os.mkdir(out_sample_dir)

            for fp in in_sample_dir.iterdir():
                out_fp = os.path.join(out_sample_dir, fp.name)
                duplicate(fp, out_fp)

    return collated_indices


def collate_alignments(alignments: BAMDirFmt) -> BAMDirFmt:
    collated_alignments = BAMDirFmt()

    for alignment in alignments:
        for _alignment in alignment.path.iterdir():
            filename = os.path.basename(_alignment)
            duplicate(_alignment, os.path.join(collated_alignments.path, filename))

    return collated_alignments


def collate_genomes(
    genomes_in: Union[DNAFASTAFormat, GenomeSequencesDirectoryFormat]
) -> GenomeSequencesDirectoryFormat:
    genomes_dir = GenomeSequencesDirectoryFormat()
    if isinstance(genomes_in[0], DNAFASTAFormat):
        for genome_file in genomes_in:
            for genome in genome_file.view(DNAIterator):
                name = genome.metadata["description"].split(",")[0]
                with open(os.path.join(genomes_dir.path, name + ".fasta"), "w") as f:
                    skbio.io.write(genome, format="fasta", into=f)
    else:
        for genome in genomes_in:
            for genome_fp in genome.path.iterdir():
                shutil.copyfile(
                    genome_fp,
                    os.path.join(genomes_dir.path, os.path.basename(genome_fp)),
                )

    return genomes_dir
