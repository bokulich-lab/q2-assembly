# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile
import unittest
from pathlib import Path
from subprocess import CalledProcessError
from unittest.mock import ANY, call, patch

from q2_types.bowtie2 import Bowtie2IndexDirFmt
from q2_types.per_sample_sequences import ContigSequencesDirFmt, MultiMAGSequencesDirFmt
from qiime2 import Artifact
from qiime2.plugin.testing import TestPluginBase
from qiime2.sdk.parallel_config import ParallelConfig

from q2_assembly.bowtie2.indexing import _index_contigs, _index_seqs, index_mags


class TestBowtie2Indexing(TestPluginBase):
    package = "q2_assembly.bowtie2.tests"

    def setUp(self):
        super().setUp()
        self.index_contigs = self.plugin.pipelines["index_contigs"]
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

    @patch("q2_assembly.bowtie2.indexing._assert_inputs_not_empty")
    @patch("subprocess.run")
    @patch("os.makedirs")
    def test_index_seqs_contigs(self, p1, p2, p3):
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

    @patch("q2_assembly.bowtie2.indexing._assert_inputs_not_empty")
    @patch("subprocess.run")
    @patch("os.makedirs")
    def test_index_seqs_mags(self, p1, p2, p3):
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

    @patch("q2_assembly.bowtie2.indexing._assert_inputs_not_empty")
    @patch(
        "subprocess.run", side_effect=CalledProcessError(returncode=123, cmd="some cmd")
    )
    @patch("os.makedirs")
    def test_index_seqs_with_error(self, p1, p2, p3):
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

    def test_index_contigs_parallel(self):
        input_contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        input_artifact = Artifact.import_data("SampleData[Contigs]", input_contigs)

        with ParallelConfig():
            (out,) = self.index_contigs.parallel(
                input_artifact,
                large_index=True,
                bmax=11,
                bmaxdivn=4,
                dcv=1024,
                offrate=5,
                ftabchars=10,
                threads=1,
            )._result()

        out.validate()
        self.assertIs(out.format, Bowtie2IndexDirFmt)

    def test_index_contigs_parallel_too_many_partitions(self):
        input_contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        input_artifact = Artifact.import_data("SampleData[Contigs]", input_contigs)
        A_MODEST_NUMBER_OF_PARTITIONS = 100000000

        with self.assertWarnsRegex(
            UserWarning, f"You have requested.*{A_MODEST_NUMBER_OF_PARTITIONS}.*2"
        ):
            with ParallelConfig():
                (out,) = self.index_contigs.parallel(
                    input_artifact,
                    large_index=True,
                    bmax=11,
                    bmaxdivn=4,
                    dcv=1024,
                    offrate=5,
                    ftabchars=10,
                    threads=1,
                    num_partitions=A_MODEST_NUMBER_OF_PARTITIONS,
                )._result()

        out.validate()
        self.assertIs(out.format, Bowtie2IndexDirFmt)

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

    def test_empty_fasta_input(self):
        with self.assertRaisesRegex(
            ValueError, r".*files were empty.*empty_contigs.fa.*second_empty_contigs.fa"
        ):
            with tempfile.TemporaryDirectory() as tempdir:
                _index_seqs(
                    fasta_fps=[
                        self.get_data_path(
                            Path("empty_contigs") / "sample1_contigs.fa"
                        ),
                        self.get_data_path(Path("empty_contigs") / "empty_contigs.fa"),
                        self.get_data_path(
                            Path("empty_contigs") / "second_empty_contigs.fa"
                        ),
                    ],
                    result_fp=Path(tempdir) / "out",
                    common_args=self.test_params_list,
                    input_type="contigs",
                )


if __name__ == "__main__":
    unittest.main()
