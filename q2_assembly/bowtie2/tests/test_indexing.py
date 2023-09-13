# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from subprocess import CalledProcessError
from unittest.mock import ANY, call, patch

from q2_types_genomics.per_sample_data import (
    ContigSequencesDirFmt,
    MultiMAGSequencesDirFmt,
)
from qiime2.plugin.testing import TestPluginBase

from q2_assembly.bowtie2.indexing import _index_contigs, _index_seqs, index_mags


class TestBowtie2Indexing(TestPluginBase):
    package = "q2_assembly.bowtie2.tests"

    def setUp(self):
        super().setUp()
        self.test_params_list = [
            "--large-index",
            "--bmax",
            "11",
            "--bmaxdivn",
            "4",
            "--dcv",
            "1024",
            "--offrate",
            "5",
            "--ftabchars",
            "10",
            "--threads",
            "1",
        ]

    @patch("subprocess.run")
    @patch("os.makedirs")
    def test_index_seqs_contigs(self, p1, p2):
        _index_seqs(
            fasta_fps=["/here/samp1_contigs.fa", "/here/samp2_contigs.fa"],
            result_fp="/there/",
            common_args=self.test_params_list,
            input_type="contigs",
        )

        p1.assert_has_calls([call("/there/samp1"), call("/there/samp2")])
        p2.assert_has_calls(
            [
                call(
                    [
                        "bowtie2-build",
                        "--large-index",
                        "--bmax",
                        "11",
                        "--bmaxdivn",
                        "4",
                        "--dcv",
                        "1024",
                        "--offrate",
                        "5",
                        "--ftabchars",
                        "10",
                        "--threads",
                        "1",
                        "/here/samp1_contigs.fa",
                        "/there/samp1/index",
                    ],
                    check=True,
                ),
                call(
                    [
                        "bowtie2-build",
                        "--large-index",
                        "--bmax",
                        "11",
                        "--bmaxdivn",
                        "4",
                        "--dcv",
                        "1024",
                        "--offrate",
                        "5",
                        "--ftabchars",
                        "10",
                        "--threads",
                        "1",
                        "/here/samp2_contigs.fa",
                        "/there/samp2/index",
                    ],
                    check=True,
                ),
            ]
        )

    @patch("subprocess.run")
    @patch("os.makedirs")
    def test_index_seqs_mags(self, p1, p2):
        _index_seqs(
            fasta_fps=["/here/smp1/mag1.fa", "/here/smp1/mag2.fa"],
            result_fp="/there/",
            common_args=self.test_params_list,
            input_type="mags",
        )

        p1.assert_has_calls([call("/there/smp1/mag1"), call("/there/smp1/mag2")])
        p2.assert_has_calls(
            [
                call(
                    [
                        "bowtie2-build",
                        "--large-index",
                        "--bmax",
                        "11",
                        "--bmaxdivn",
                        "4",
                        "--dcv",
                        "1024",
                        "--offrate",
                        "5",
                        "--ftabchars",
                        "10",
                        "--threads",
                        "1",
                        "/here/smp1/mag1.fa",
                        "/there/smp1/mag1/index",
                    ],
                    check=True,
                ),
                call(
                    [
                        "bowtie2-build",
                        "--large-index",
                        "--bmax",
                        "11",
                        "--bmaxdivn",
                        "4",
                        "--dcv",
                        "1024",
                        "--offrate",
                        "5",
                        "--ftabchars",
                        "10",
                        "--threads",
                        "1",
                        "/here/smp1/mag2.fa",
                        "/there/smp1/mag2/index",
                    ],
                    check=True,
                ),
            ]
        )

    @patch(
        "subprocess.run", side_effect=CalledProcessError(returncode=123, cmd="some cmd")
    )
    @patch("os.makedirs")
    def test_index_seqs_with_error(self, p1, p2):
        with self.assertRaisesRegex(
            Exception, "An error.*while running Bowtie2.*code 123"
        ):
            _index_seqs(
                fasta_fps=["/here/samp1/mag1.fa"],
                result_fp="/there/",
                common_args=self.test_params_list,
                input_type="mags",
            )

    @patch("q2_assembly.bowtie2.indexing._index_seqs")
    def test_index_contigs(self, p):
        input_contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        _index_contigs(
            input_contigs,
            large_index=True,
            bmax=11,
            bmaxdivn=4,
            dcv=1024,
            offrate=5,
            ftabchars=10,
            threads=1,
        )

        exp_contigs = [f"{str(input_contigs)}/sample{x+1}_contigs.fa" for x in range(2)]
        p.assert_called_with(exp_contigs, ANY, self.test_params_list, "contigs")

    @patch("q2_assembly.bowtie2.indexing._index_seqs")
    def test_index_mags(self, p):
        input_mags = MultiMAGSequencesDirFmt(self.get_data_path("mags"), "r")
        index_mags(
            input_mags,
            large_index=True,
            bmax=11,
            bmaxdivn=4,
            dcv=1024,
            offrate=5,
            ftabchars=10,
            threads=1,
        )

        exp_mags = [
            f"{str(input_mags)}/sample{x+1}/mag{y+1}.fa"
            for x in range(2)
            for y in range(2)
        ]
        p.assert_called_with(exp_mags, ANY, self.test_params_list, "mags")


if __name__ == "__main__":
    unittest.main()
