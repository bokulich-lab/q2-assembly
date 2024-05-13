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
OPTIONAL_COLS_MAP = {
    "# contigs (>= 1000 bp)": "no_contigs_1000",
    "# contigs (>= 5000 bp)": "no_contigs_5000",
    "# contigs (>= 10000 bp)": "no_contigs_10000",
    "# contigs (>= 25000 bp)": "no_contigs_25000",
    "# contigs (>= 50000 bp)": "no_contigs_50000",
    "Total length (>= 0 bp)": "total_length_0",
    "Total length (>= 1000 bp)": "total_length_1000",
    "Total length (>= 5000 bp)": "total_length_5000",
    "Total length (>= 10000 bp)": "total_length_10000",
    "Total length (>= 25000 bp)": "total_length_25000",
    "Total length (>= 50000 bp)": "total_length_50000",
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

RESHUFFLED_COLUMNS = [
    "id",
    "total_length",
    "no_contigs_0",
    "no_contigs",
    "longest_contig",
    "n50",
    "l50",
    "n90",
    "l90",
]
