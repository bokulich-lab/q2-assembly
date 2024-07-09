# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

MANDATORY_COLS_MAP = {
    "Assembly": "id",
    "# contigs (>= 0 bp)": "no_contigs_0",
    "# contigs": "no_contigs",
    "Largest contig": "longest_contig",
    "Total length": "total_length",
    "N50": "n50",
    "N90": "n90",
    "L50": "l50",
    "L90": "l90",
}


def initialize_optional_cols_map(contig_thresholds: list):
    """
    This function will initialize the optional_cols_map dictionary
    with the contig thresholds provided by the user.

    Args:
        - contig_thresholds(list): list of contig thresholds

    Returns:
        optional_cols_map(dict): dictionary with the optional columns
    """
    optional_cols_map = {
        "Total length (>= 0 bp)": "total_length_0",
        "Reference length": "reference_length",
        "auN": "aun",
        "# misassemblies": "no_misassemblies",
        "# misassembled contigs": "no_misassembled_contigs",
        "Misassembled contigs length": "misassembled_contigs_length",
        "# local misassemblies": "no_local_misassemblies",
        "# scaffold gap ext. mis.": "no_scaffold_gap_ext_mis",
        "# scaffold gap loc. mis.": "no_scaffold_gap_loc_mis",
        "# unaligned mis. contigs": "no_unaligned_mis_contigs",
        "# unaligned contigs": "no_unaligned_contigs",
        "Unaligned length": "unaligned_length",
        "Genome fraction (%)": "genome_fraction_percentage",
        "Duplication ratio": "duplication_ratio",
        "# N's per 100 kbp": "no_ns_per100_kbp",
        "# mismatches per 100 kbp": "no_mismatches_per100_kbp",
        "# indels per 100 kbp": "no_indels_per100_kbp",
        "Largest alignment": "largest_alignment",
        "Total aligned length": "total_aligned_length",
        "auNA": "auna",
        "NA50": "na50",
        "NA90": "na90",
        "LA50": "la50",
        "LA90": "la90",
    }
    for threshold in contig_thresholds:
        optional_cols_map[
            f"Total length (>= {threshold} bp)"
        ] = f"total_length_{threshold}"
        if threshold != 0:  # this is included in mandatory cols
            optional_cols_map[
                f"# contigs (>= {threshold} bp)"
            ] = f"no_contigs_{threshold}"

    return optional_cols_map
