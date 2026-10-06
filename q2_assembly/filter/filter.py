# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import warnings

import shortuuid
import skbio
from q2_types.feature_data import FeatureData
from q2_types.feature_data_mag import Contig
from q2_types.per_sample_sequences import ContigSequencesDirFmt
from qiime2 import Metadata
from qiime2.plugin import Properties
from qiime2.util import duplicate


def _get_metadata_ids(metadata: Metadata = None, where: str = None) -> set:
    if metadata is None:
        return set()
    metadata_ids = metadata.get_ids(where=where)
    if not metadata_ids:
        print("The filter query returned no IDs to filter out.")
    return metadata_ids


def _check_missing_ids(ids: list, present_ids: set):
    missing_ids = set(ids or []) - present_ids
    if missing_ids:
        raise ValueError(
            "The following IDs are not present in the contig data: "
            f"{', '.join(sorted(missing_ids))}"
        )


# Renamed contig IDs are assumed to be '<sample_id><separator><shortuuid>',
SHORTUUID_LENGTH = shortuuid.ShortUUID().encoded_length()


def _match_contig_ids(contig_id: str, selected_ids: set[str]) -> set[str]:
    candidates = {contig_id, contig_id[-SHORTUUID_LENGTH:]}
    return candidates & selected_ids


def _match_sample_prefix_ids(contig_id: str, selected_ids: set[str]) -> set[str]:
    suffix_length = 1 + SHORTUUID_LENGTH
    if len(contig_id) < suffix_length:
        return set()
    candidates = {contig_id[:-suffix_length]}
    return candidates & selected_ids


def _filter_fasta(
    in_fp: str,
    out_fp: str,
    on: str,
    length_threshold: int,
    selected_ids: set,
    exclude_ids: bool,
) -> tuple:

    match_ids = (
        _match_sample_prefix_ids if on == "sample_prefix" else _match_contig_ids
    )
    kept, removed, found_ids = 0, 0, set()
    with open(out_fp, "w") as f_out:
        for contig in skbio.io.read(in_fp, format="fasta"):
            # do not filter contigs by ID, length threshold only applies
            if selected_ids is None:
                contig_selected = True
            else:
                _ids = match_ids(contig.metadata["id"], selected_ids)
                found_ids.update(_ids)
                if exclude_ids:
                    contig_selected = not bool(_ids)
                else:
                    contig_selected = bool(_ids)

            if contig_selected and len(contig) >= length_threshold:
                skbio.io.write(contig, format="fasta", into=f_out)
                kept += 1
            else:
                removed += 1
    return kept, removed, found_ids


def _write_filtered_contig_file(
    file_id: str,
    in_fp: str,
    out_fp: str,
    on: str,
    length_threshold: int,
    selected_ids: set[str] | None,
    exclude_ids: bool,
) -> set[str]:

    if on == "sample" and length_threshold == 0:
        duplicate(in_fp, out_fp)
        return set()

    # on-sample remaining part is only length filtering
    contig_ids = None if on == "sample" else selected_ids
    kept, removed, found_ids = _filter_fasta(
        in_fp, out_fp, on, length_threshold, contig_ids, exclude_ids
    )
    print(
        f"File {file_id}: {removed + kept} contigs\n  {removed} "
        f"contigs removed\n  {kept} contigs retained"
    )
    return found_ids


def _validate_filtered_contigs(contigs: ContigSequencesDirFmt, on: str):
    contig_files = contigs.sample_dict()
    # only reachable for SampleData[Contig] branch
    if not contig_files:
        if on == "sample":
            raise ValueError("No samples remain after filtering.")
        raise ValueError("No contigs remain after filtering.")

    if all(os.path.getsize(fp) == 0 for fp in contig_files.values()):
        warnings.warn(
            "No contigs remain after filtering - the output contains only "
            "empty files.",
            UserWarning,
        )


def _filter_contigs(
    contigs: ContigSequencesDirFmt,
    on: str,
    ids: list = None,
    metadata: Metadata = None,
    where: str = None,
    exclude_ids: bool = False,
    length_threshold: int = 0,
    remove_empty: bool = False,
) -> ContigSequencesDirFmt:
    metadata_ids = _get_metadata_ids(metadata, where)

    contig_files = contigs.sample_dict()

    if on == "sample":
        samples_to_keep = set(contig_files)
        _check_missing_ids(ids, samples_to_keep)

    selected_ids = None
    if ids is not None or metadata is not None:
        selected_ids = metadata_ids.union(ids or [])

    if on == "sample" and selected_ids is not None:
        if exclude_ids:
            samples_to_keep -= selected_ids
        else:
            samples_to_keep &= selected_ids

    if length_threshold > 0:
        print(
            f"Filtering contigs by length - only contigs >= {length_threshold} "
            f"bp long will be retained."
        )

    results = ContigSequencesDirFmt()
    empty_files = []
    found_ids = set()
    for file_id, file_fp in contig_files.items():
        if on == "sample" and file_id not in samples_to_keep:
            continue

        out_fp = os.path.join(str(results), os.path.basename(file_fp))
        file_found_ids = _write_filtered_contig_file(
            file_id,
            file_fp,
            out_fp,
            on=on,
            length_threshold=length_threshold,
            selected_ids=selected_ids,
            exclude_ids=exclude_ids,
        )
        found_ids.update(file_found_ids)
        if remove_empty and os.path.getsize(out_fp) == 0:
            os.remove(out_fp)
            empty_files.append(file_id)

    # Missing IDs can only be determined after reading all input FASTA files.
    if on != "sample":
        _check_missing_ids(ids, found_ids)

    if empty_files:
        print(f"Removing empty samples: {', '.join(sorted(empty_files))}")
    _validate_filtered_contigs(results, on)

    return results


def filter_contigs(
    ctx,
    contigs,
    on="sample",
    ids=None,
    metadata=None,
    where=None,
    exclude_ids=False,
    length_threshold=0,
    remove_empty=False,
):
    if not any([ids, metadata, length_threshold, remove_empty]):
        raise ValueError(
            "At least one of the following parameters must be provided: "
            "ids, metadata, length_threshold, remove_empty."
        )

    if where is not None and metadata is None:
        raise ValueError("Metadata must be provided if 'where' is specified.")

    if exclude_ids and ids is None and metadata is None:
        raise ValueError(
            "Either 'ids' or metadata must be provided if 'exclude_ids' is True."
        )

    if (
        contigs.type <= FeatureData[Contig % Properties("pooled")]
        and on == "sample"
    ):
        # sample IDs are matched against the contig ID prefix
        on = "sample_prefix"

    _filter_contigs = ctx.get_action("assembly", "_filter_contigs")
    (filtered_contigs,) = _filter_contigs(
        contigs=contigs,
        on=on,
        ids=ids,
        metadata=metadata,
        where=where,
        exclude_ids=exclude_ids,
        length_threshold=length_threshold,
        remove_empty=remove_empty,
    )
    return filtered_contigs
