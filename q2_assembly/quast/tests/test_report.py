# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import contextlib
import tempfile
import unittest

from qiime2.plugin.testing import TestPluginBase

from q2_assembly.quast.report import initialize_optional_cols_map

actual_optional_cols_map_1 = {
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
actual_optional_cols_map_2 = {
    "# contigs (>= 1000 bp)": "no_contigs_1000",
    "# contigs (>= 5000 bp)": "no_contigs_5000",
    "# contigs (>= 100000 bp)": "no_contigs_100000",
    "# contigs (>= 250000 bp)": "no_contigs_250000",
    "# contigs (>= 500000 bp)": "no_contigs_500000",
    "Total length (>= 0 bp)": "total_length_0",
    "Total length (>= 1000 bp)": "total_length_1000",
    "Total length (>= 5000 bp)": "total_length_5000",
    "Total length (>= 100000 bp)": "total_length_100000",
    "Total length (>= 250000 bp)": "total_length_250000",
    "Total length (>= 500000 bp)": "total_length_500000",
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


class TestQuastUtils(TestPluginBase):
    package = "q2_assembly.tests"

    def setUp(self):
        super().setUp()
        with contextlib.ExitStack() as stack:
            self._tmp = stack.enter_context(tempfile.TemporaryDirectory())
            self.addCleanup(stack.pop_all().close)

    def test_initialize_optional_cols_map(self):
        contig_thresholds = [1000, 5000, 10000, 25000, 50000]
        created_optional_cols_map = initialize_optional_cols_map(contig_thresholds)

        assert set(created_optional_cols_map.keys()) == set(
            actual_optional_cols_map_1.keys()
        )
        assert set(created_optional_cols_map.values()) == set(
            actual_optional_cols_map_1.values()
        )

    def test_intialize_optional_cols_map_2(self):
        contig_thresholds = [1000, 5000, 100000, 250000, 500000]
        created_optional_cols_map = initialize_optional_cols_map(contig_thresholds)

        assert set(created_optional_cols_map.keys()) == set(
            actual_optional_cols_map_2.keys()
        )
        assert set(created_optional_cols_map.values()) == set(
            actual_optional_cols_map_2.values()
        )


if __name__ == "__main__":
    unittest.main()
