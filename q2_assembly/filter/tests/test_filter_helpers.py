# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from unittest import TestCase

from q2_assembly.filter.filter import _match_contig_ids, _match_sample_prefix_ids


class TestContigIdMatching(TestCase):
    def test_sample_prefix_separator_in_sample_id(self):
        self.assertSetEqual(
            _match_sample_prefix_ids(
                "my_sample_1_mtjebimcR24S9DZ62TY6Fh", {"my_sample_1"}
            ),
            {"my_sample_1"},
        )

    def test_sample_prefix_not_matched_by_shorter_sample_id(self):
        self.assertSetEqual(
            _match_sample_prefix_ids("sample10:mtjebimcR24S9DZ62TY6Fh", {"sample1"}),
            set(),
        )

    def test_sample_prefix_unrecognized_contig_id(self):
        self.assertSetEqual(
            _match_sample_prefix_ids("k141_19", {"sample1", "k141"}), set()
        )

    def test_contig_not_matched_by_partial_suffix(self):
        self.assertSetEqual(_match_contig_ids("k141_19", {"9", "k141_9"}), set())

    def test_contig_matches_both_whole_id_and_shortuuid(self):
        contig_id = "sample1:mtjebimcR24S9DZ62TY6Fh"
        selected_ids = {contig_id, "mtjebimcR24S9DZ62TY6Fh", "sample1", "unrelated"}
        self.assertSetEqual(
            _match_contig_ids(contig_id, selected_ids),
            {contig_id, "mtjebimcR24S9DZ62TY6Fh"},
        )
