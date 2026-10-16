# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from unittest import TestCase

from parameterized import parameterized

from q2_assembly._action_params import ALLOWED_SEPARATORS
from q2_assembly.filter.filter import _match_contig_ids, _match_sample_prefix_ids


class TestContigIdParsing(TestCase):

    @parameterized.expand(
        [
            ("shortuuid", "mtjebimcR24S9DZ62TY6Fh"),
            ("uuid3", "2f1f2fae-d96e-3481-adab-21285b467f10"),
            ("uuid4", "550e8400-e29b-41d4-a716-446655440000"),
            ("uuid5", "fbf870ed-9777-513b-851c-3bbd7b3a10f8"),
        ]
    )
    def test_sample_prefix_any_separator(self, uuid_type, contig_uuid):
        sample_id = "my_sample:C.1"
        for separator in ALLOWED_SEPARATORS:
            with self.subTest(separator=separator):
                contig_id = f"{sample_id}{separator}{contig_uuid}"
                self.assertSetEqual(
                    _match_sample_prefix_ids(contig_id, {sample_id, "my_sample"}),
                    {sample_id},
                )

    @parameterized.expand(
        [
            ("shortuuid", "mtjebimcR24S9DZ62TY6Fh"),
            ("uuid3", "2f1f2fae-d96e-3481-adab-21285b467f10"),
            ("uuid4", "550e8400-e29b-41d4-a716-446655440000"),
            ("uuid5", "fbf870ed-9777-513b-851c-3bbd7b3a10f8"),
        ]
    )
    def test_contig_any_separator(self, uuid_type, contig_uuid):
        for separator in ALLOWED_SEPARATORS:
            with self.subTest(separator=separator):
                contig_id = f"my_sample:C.1{separator}{contig_uuid}"
                self.assertSetEqual(
                    _match_contig_ids(
                        contig_id, {contig_id, contig_uuid, contig_uuid[-10:]}
                    ),
                    {contig_id, contig_uuid},
                )

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
