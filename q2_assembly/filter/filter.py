# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os

import skbio
from q2_types.per_sample_sequences import ContigSequencesDirFmt
from qiime2 import Metadata
from qiime2.util import duplicate


def _find_empty_samples(samples: dict) -> set:
    empty_samples = set()
    for sample_id, sample_fp in samples.items():
        if os.path.getsize(sample_fp) == 0:
            empty_samples.add(sample_id)
    return empty_samples


def _filter_by_length(
    contigs: ContigSequencesDirFmt, threshold: int
) -> ContigSequencesDirFmt:
    results = ContigSequencesDirFmt()
    print(
        f"Filtering contigs by length - only contigs >= {threshold} bp long will "
        f"be retained."
    )
    for sample_id, sample_fp in contigs.sample_dict().items():
        out_fp = os.path.join(str(results), f"{sample_id}_contigs.fa")
        keep, remove = 0, 0
        with open(out_fp, "w") as f_out:
            for contig in skbio.io.read(sample_fp, format="fasta"):
                if len(contig) >= threshold:
                    skbio.io.write(contig, format="fasta", into=f_out)
                    keep += 1
                else:
                    remove += 1
            print(
                f"Sample {sample_id}: {remove + keep} contigs\n  {remove} contigs "
                f"removed\n  {keep} contigs retained"
            )
    return results


def filter_contigs(
    contigs: ContigSequencesDirFmt,
    metadata: Metadata = None,
    where: str = None,
    exclude_ids: bool = False,
    length_threshold: int = 0,
    remove_empty: bool = False,
) -> ContigSequencesDirFmt:
    if not any([metadata, length_threshold, remove_empty]):
        raise ValueError(
            "At least one of the following parameters must be provided: "
            "metadata, length_threshold, remove_empty."
        )

    if metadata and not where:
        raise ValueError(
            "A filter query must be provided through the 'when' parameter "
            "when filtering by metadata."
        )

    if length_threshold > 0:
        contigs = _filter_by_length(contigs, length_threshold)

    results = ContigSequencesDirFmt()
    samples = contigs.sample_dict()
    ids_to_keep = set(samples.keys())
    if remove_empty:
        ids_to_remove = _find_empty_samples(samples)
        ids_to_keep -= ids_to_remove
        if ids_to_remove:
            print(f"Removing empty samples: {', '.join(sorted(ids_to_remove))}")

    if metadata:
        selected_ids = metadata.get_ids(where=where)
        if not selected_ids:
            print("The filter query returned no IDs to filter out.")

        if exclude_ids:
            ids_to_keep -= set(selected_ids)
        else:
            ids_to_keep &= set(selected_ids)

    if len(ids_to_keep) == 0:
        raise ValueError("No samples remain after filtering.")

    try:
        for _id in ids_to_keep:
            duplicate(samples[_id], os.path.join(str(results), f"{_id}_contigs.fa"))
    except KeyError:
        raise ValueError(f"{_id!r} is not a sample present in the contig data.")

    return results
