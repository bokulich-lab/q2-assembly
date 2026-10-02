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
from q2_types.per_sample_sequences import ContigSequencesDirFmt, Contigs
from q2_types.sample_data import SampleData
from qiime2 import Metadata
from qiime2.plugin import Properties
from qiime2.util import duplicate


def _get_metadata_ids(metadata: Metadata = None, where: str = None) -> set:
    if metadata is None:
        return set()
    metadata_ids = set(metadata.get_ids(where=where))
    if not metadata_ids:
        print("The filter query returned no IDs to filter out.")
    return metadata_ids


def _check_missing_ids(ids: list, metadata_ids: set, present_ids: set):
    missing_ids = set(ids or []) - present_ids
    if missing_ids:
        raise ValueError(
            "The following IDs are not present in the contig data: "
            f"{', '.join(sorted(missing_ids))}"
        )
    if metadata_ids and not metadata_ids & present_ids:
        raise ValueError(
            "None of the metadata IDs are present in the contig data. Check "
            "that `on` parameter matches the metadata ID type (sample or contig)."
        )


# Renamed contig IDs are assumed to be '<sample_id><separator><shortuuid>',
# Contigs renamed with uuid3/uuid4/uuid5 are not matched correctly.
SHORTUUID_LENGTH = shortuuid.ShortUUID().encoded_length()


def _matched_ids(contig_id: str, on: str, selected_ids: set) -> set:

    if on == "sample_prefix":
        return {
            _id
            for _id in selected_ids
            if contig_id.startswith(_id)
            and len(contig_id) == len(_id) + 1 + SHORTUUID_LENGTH
        }
    return {
        _id
        for _id in selected_ids
        if contig_id.endswith(_id) and len(_id) in (len(contig_id), SHORTUUID_LENGTH)
    }


def _is_selected(_ids: set, selected_ids: set, exclude_ids: bool) -> bool:
    return selected_ids is None or bool(_ids & selected_ids) != exclude_ids


def _filter_fasta(
    in_fp: str,
    out_fp: str,
    on: str,
    length_threshold: int,
    selected_ids: set,
    exclude_ids: bool,
) -> tuple:

    kept, removed, found_ids = 0, 0, set()
    with open(out_fp, "w") as f_out:
        for contig in skbio.io.read(in_fp, format="fasta"):
            # no selection by contig IDs, only the length filter applies
            if selected_ids is None:
                contig_selected = True
            else:
                _ids = _matched_ids(contig.metadata["id"], on, selected_ids)
                found_ids.update(_ids)
                contig_selected = _is_selected(_ids, selected_ids, exclude_ids)

            if contig_selected and len(contig) >= length_threshold:
                skbio.io.write(contig, format="fasta", into=f_out)
                kept += 1
            else:
                removed += 1
    return kept, removed, found_ids


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
        _check_missing_ids(ids, metadata_ids, set(contig_files))

    selected_ids = None
    if ids is not None or metadata is not None:
        selected_ids = metadata_ids.union(ids or [])

    if length_threshold > 0:
        print(
            f"Filtering contigs by length - only contigs >= {length_threshold} "
            f"bp long will be retained."
        )

    results = ContigSequencesDirFmt()
    empty_files = []
    found_ids = set()
    for file_id, file_fp in contig_files.items():
        # non-selected samples are skipped
        if on == "sample" and not _is_selected(
            {file_id}, selected_ids, exclude_ids
        ):
            continue

        out_fp = os.path.join(str(results), f"{file_id}.fa")
        # the only case where contig files are ignored       
        if on == "sample" and length_threshold == 0:
            is_empty = os.path.getsize(file_fp) == 0
            # empty samples to be removed are not copied at all
            if remove_empty and is_empty:
                empty_files.append(file_id)
                continue
            duplicate(file_fp, out_fp)
        # filtering applied inside contig files > _filter_fasta
        else:
            # on="sample": samples were already selected above, so only the
            # length filter applies within the file
            contig_ids = None if on == "sample" else selected_ids
            kept, removed, file_found_ids = _filter_fasta(
                file_fp, out_fp, on, length_threshold, contig_ids, exclude_ids
            )
            found_ids.update(file_found_ids)
            print(
                f"Sample {file_id}: {removed + kept} contigs\n  {removed} "
                f"contigs removed\n  {kept} contigs retained"
            )
            # emptiness is only known after filtering, so the file is removed
            if remove_empty and kept == 0:
                os.remove(out_fp)
                empty_files.append(file_id)

    if on != "sample":
        _check_missing_ids(ids, metadata_ids, found_ids)

    # samples excluded or removed as empty - only reachable for SampleData
    if empty_files:
        print(f"Removing empty samples: {', '.join(sorted(empty_files))}")
    if not results.sample_dict():
        if on == "sample":
            raise ValueError("No samples remain after filtering.")
        raise ValueError("No contigs remain after filtering.")

    # empty files were kept (remove_empty=False)
    if all(os.path.getsize(fp) == 0 for fp in results.sample_dict().values()):
        warnings.warn(
            "No contigs remain after filtering - the output contains only "
            "empty files.",
            UserWarning,
        )

    return results


def filter_contigs(
    ctx,
    contigs,
    on="None",
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

    if contigs.type <= SampleData[Contigs]:
        on = "sample" if on == "None" else on

    elif contigs.type <= FeatureData[Contig % Properties("pooled")]:
        # sample IDs are matched against the contig ID prefix
        on = "sample_prefix" if on == "sample" else "contig"

    else:  # FeatureData[Contig % Properties("coassembly")]
        # `on="sample"` is not allowed by the TypeMap
        on = "contig"

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
