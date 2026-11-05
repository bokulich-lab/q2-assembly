# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import shutil

from q2_types.bowtie2 import Bowtie2IndexDirFmt
from q2_types.per_sample_sequences import BAMDirFmt, ContigSequencesDirFmt
from qiime2.util import duplicate

from q2_assembly._utils import _split_contig_id, modify_contig_ids, run_command


def rename_contigs(
    contigs: ContigSequencesDirFmt,
    uuid_type: str,
    include_sample_id: bool = True,
    separator: str = ":",
) -> ContigSequencesDirFmt:
    renamed_contigs = ContigSequencesDirFmt()
    sep = separator if include_sample_id else None

    for sample_id, contigs_fp in contigs.sample_dict().items():
        out_fp = os.path.join(renamed_contigs.path, os.path.basename(contigs_fp))
        shutil.copyfile(contigs_fp, out_fp)
        modify_contig_ids(out_fp, sample_id, uuid_type, sep)

    return renamed_contigs


def _check_contig_id(contig_id: str, sample_id: str):
    parsed_id = _split_contig_id(contig_id)
    if parsed_id is None:
        raise ValueError(
            f'Contig ID "{contig_id}" in sample "{sample_id}" is not of the form '
            f'"{sample_id}<separator><uuid>". Please run rename-contigs with '
            "include_sample_id enabled before pooling."
        )
    if parsed_id[0] != sample_id:
        raise ValueError(
            f'Contig ID "{contig_id}" in sample "{sample_id}" is prefixed with '
            f'a different sample ID ("{parsed_id[0]}").'
        )


def pool_contigs(contigs: ContigSequencesDirFmt) -> ContigSequencesDirFmt:
    """Join contigs from one or more samples into one pooled set."""
    pooled_contigs = ContigSequencesDirFmt()
    out_fp = os.path.join(pooled_contigs.path, "pooled.fa")

    with open(out_fp, "w") as f_out:
        for sample_id, contigs_fp in contigs.sample_dict().items():
            with open(contigs_fp) as f_in:
                for line in f_in:
                    if line.startswith(">"):
                        contig_id = (line[1:].split() or [""])[0]
                        _check_contig_id(contig_id, sample_id)
                    f_out.write(line if line.endswith("\n") else line + "\n")

    return pooled_contigs


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


def collate_alignments(alignment_maps: BAMDirFmt) -> BAMDirFmt:
    collated_alignments = BAMDirFmt()

    for alignment in alignment_maps:
        for _alignment in alignment.path.iterdir():
            filename = os.path.basename(_alignment)
            duplicate(_alignment, os.path.join(collated_alignments.path, filename))

    return collated_alignments


def sort_alignment_maps(
    alignment_maps: BAMDirFmt,
) -> BAMDirFmt:
    file_dict = alignment_maps.file_dict(relative=True)
    out_dir = BAMDirFmt()

    for fp in file_dict.values():
        run_command(
            [
                "samtools",
                "sort",
                str(alignment_maps.path / fp),
                "-o",
                str(out_dir.path / fp),
            ],
            verbose=True,
        )

    return out_dir
