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
from pathlib import Path

from q2_types.bowtie2 import Bowtie2IndexDirFmt
from q2_types.per_sample_sequences import BAMDirFmt, ContigSequencesDirFmt
from qiime2.util import duplicate

from q2_assembly._utils import modify_contig_ids, run_command


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
    map_fps = glob.glob(os.path.join(str(alignment_maps), "*.bam"))
    out_dir = BAMDirFmt()

    for fp in map_fps:
        samp_name = Path(fp).stem.rsplit("_alignment", 1)[0]
        sorted_bam = os.path.join(str(out_dir), f"{samp_name}_alignment_sorted.bam")
        run_command(["samtools", "sort", str(fp), "-o", sorted_bam], verbose=True)

    return out_dir
