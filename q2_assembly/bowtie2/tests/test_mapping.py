# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import unittest
from subprocess import CalledProcessError
from unittest.mock import ANY, call, patch

import pandas as pd
from q2_types.bowtie2 import Bowtie2IndexDirFmt
from q2_types.per_sample_sequences import (
    BAMDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from qiime2 import Artifact
from qiime2.plugin.testing import TestPluginBase
from qiime2.sdk.parallel_config import ParallelConfig

from q2_assembly._utils import get_relative_data_path
from q2_assembly.bowtie2.mapping import (
    _gather_sample_data,
    _map_reads_to_contigs,
    _map_sample_reads,
)


class MockTempDir(tempfile.TemporaryDirectory):
    pass


class TestBowtie2Mapping(TestPluginBase):
    package = "q2_assembly.bowtie2.tests"
    root_test_package = "q2_assembly.tests"

    def setUp(self):
        super().setUp()
        self.map_reads_to_contigs = self.plugin.pipelines["map_reads_to_contigs"]
        self.test_params_list = [
            "--trim5",
            "10",
            "-N",
            "1",
            "-L",
            "22",
            "-i",
            "L,1,0.5",
            "--n-ceil",
            "L,0,0.15",
            "--dpad",
            "15",
            "--gbar",
            "4",
            "--ma",
            "2",
            "--mp",
            "6",
            "--np",
            "1",
            "--rdg",
            "5,3",
            "--rfg",
            "5,3",
            "-D",
            "15",
            "-R",
            "2",
            "--maxins",
            "500",
            "--ff",
            "--threads",
            "1",
            "--very-fast",
        ]
        self.test_index = Bowtie2IndexDirFmt()
        self.test_result = BAMDirFmt()
        self.test_samples = {
            "sample1": {
                "fwd": "sample1_00_L001_R1_001.fastq.gz",
                "rev": "sample1_00_L001_R2_001.fastq.gz",
                "index": os.path.join(str(self.test_index), "sample1", "index"),
            },
            "sample2": {
                "fwd": "sample2_00_L001_R1_001.fastq.gz",
                "rev": "sample2_00_L001_R2_001.fastq.gz",
                "index": os.path.join(str(self.test_index), "sample2", "index"),
            },
        }

    def read_manifest_file(self, manifest_type):
        fp = self.get_data_path(f"manifest-{manifest_type}.csv")
        return pd.read_csv(fp, sep=",", index_col=0, header=0)

    def test_gather_sample_data_paired(self):
        manifest = self.read_manifest_file("paired")
        for s in ["sample1", "sample2"]:
            os.makedirs(os.path.join(str(self.test_index), s))

        obs = _gather_sample_data(
            indexed_contigs=self.test_index, reads_manifest=manifest, paired=True
        )
        exp = self.test_samples
        self.assertDictEqual(obs, exp)

    def test_gather_sample_data_single(self):
        manifest = self.read_manifest_file("single")
        for s in ["sample1", "sample2"]:
            os.makedirs(os.path.join(str(self.test_index), s))

        obs = _gather_sample_data(
            indexed_contigs=self.test_index, reads_manifest=manifest, paired=False
        )
        exp = self.test_samples
        for s in exp:
            exp[s]["rev"] = None
        self.assertDictEqual(obs, exp)

    def test_gather_sample_data_paired_missing_sample(self):
        manifest = self.read_manifest_file("wrong-samples")
        for s in ["sample1", "sample2"]:
            os.makedirs(os.path.join(str(self.test_index), s))

        with self.assertRaisesRegex(
            Exception, "Index files missing for sample sample3."
        ):
            _gather_sample_data(
                indexed_contigs=self.test_index, reads_manifest=manifest, paired=True
            )

    @patch("shutil.move")
    @patch("subprocess.run")
    @patch("tempfile.TemporaryDirectory")
    def test_map_sample_reads_paired(self, p1, p2, p3):
        test_temp_dir = MockTempDir()
        p1.return_value = test_temp_dir

        exp_bam = os.path.join(test_temp_dir.name, "alignment.bam")
        for s, s_props in self.test_samples.items():
            _map_sample_reads(
                common_args=self.test_params_list,
                paired=True,
                sample_name=s,
                sample_inputs=s_props,
                result_fp=str(self.test_result),
            )

            exp_calls = [
                call(
                    ["bowtie2"]
                    + self.test_params_list
                    + [
                        "-x",
                        os.path.join(str(self.test_index), s, "index"),
                        "-1",
                        f"{s}_00_L001_R1_001.fastq.gz",
                        "-2",
                        f"{s}_00_L001_R2_001.fastq.gz",
                    ],
                    check=True,
                    capture_output=True,
                ),
                call(["samtools", "view", "-bS", "-o", exp_bam], input=ANY, check=True),
            ]

            p2.assert_has_calls(exp_calls)
            p3.assert_called_with(
                exp_bam, os.path.join(str(self.test_result), f"{s}_alignment.bam")
            )

    @patch(
        "subprocess.run", side_effect=CalledProcessError(returncode=123, cmd="some cmd")
    )
    def test_map_sample_reads_error(self, p):
        for s, s_props in self.test_samples.items():
            with self.assertRaisesRegex(
                Exception, "An error.*while running Bowtie2.*code 123"
            ):
                _map_sample_reads(
                    common_args=self.test_params_list,
                    paired=True,
                    sample_name=s,
                    sample_inputs=s_props,
                    result_fp=str(self.test_result),
                )

    @patch("shutil.move")
    @patch("subprocess.run")
    @patch("tempfile.TemporaryDirectory")
    def test_map_sample_reads_single(self, p1, p2, p3):
        test_temp_dir = MockTempDir()
        p1.return_value = test_temp_dir

        exp_bam = os.path.join(test_temp_dir.name, "alignment.bam")
        for s, s_props in self.test_samples.items():
            s_props["rev"] = None
            _map_sample_reads(
                common_args=self.test_params_list,
                paired=False,
                sample_name=s,
                sample_inputs=s_props,
                result_fp=str(self.test_result),
            )

            exp_calls = [
                call(
                    ["bowtie2"]
                    + self.test_params_list
                    + [
                        "-x",
                        os.path.join(str(self.test_index), s, "index"),
                        "-U",
                        f"{s}_00_L001_R1_001.fastq.gz",
                    ],
                    check=True,
                    capture_output=True,
                ),
                call(["samtools", "view", "-bS", "-o", exp_bam], input=ANY, check=True),
            ]

            p2.assert_has_calls(exp_calls)
            p3.assert_called_with(
                exp_bam, os.path.join(str(self.test_result), f"{s}_alignment.bam")
            )

    @patch("q2_assembly.bowtie2.mapping._map_sample_reads")
    @patch("q2_assembly.bowtie2.mapping._gather_sample_data")
    def test_map_reads_to_contigs_paired(self, p1, p2):
        input_reads = self.get_data_path("reads/paired-end")
        input_index = self.get_data_path("indices/from_contigs")
        reads = SingleLanePerSamplePairedEndFastqDirFmt(input_reads, mode="r")
        index = Bowtie2IndexDirFmt(input_index, mode="r")

        p1.return_value = self.test_samples

        _map_reads_to_contigs(
            indexed_contigs=index,
            reads=reads,
            trim5=10,
            n=1,
            i="L,1,0.5",
            valid_mate_orientations="ff",
            mode="global",
            sensitivity="very-fast",
        )

        p1.assert_called_with(index, ANY, True)
        pd.testing.assert_frame_equal(
            p1.call_args[0][1], reads.manifest.view(pd.DataFrame)
        )

        exp_calls = []
        for s, s_props in self.test_samples.items():
            exp_calls.append(call(self.test_params_list, True, s, s_props, ANY))
        p2.assert_has_calls(exp_calls)

    @patch("q2_assembly.bowtie2.mapping._map_sample_reads")
    @patch("q2_assembly.bowtie2.mapping._gather_sample_data")
    def test_map_reads_to_contigs_single(self, p1, p2):
        input_reads = self.get_data_path("reads/single-end")
        input_index = self.get_data_path("indices/from_contigs")
        reads = SingleLanePerSampleSingleEndFastqDirFmt(input_reads, mode="r")
        index = Bowtie2IndexDirFmt(input_index, mode="r")

        for s in self.test_samples:
            self.test_samples[s]["rev"] = None
        p1.return_value = self.test_samples

        _map_reads_to_contigs(
            indexed_contigs=index,
            reads=reads,
            trim5=10,
            n=1,
            i="L,1,0.5",
            valid_mate_orientations="ff",
            mode="global",
            sensitivity="very-fast",
        )

        p1.assert_called_with(index, ANY, False)
        pd.testing.assert_frame_equal(
            p1.call_args[0][1], reads.manifest.view(pd.DataFrame)
        )

        exp_calls = []
        for s, s_props in self.test_samples.items():
            exp_calls.append(call(self.test_params_list, False, s, s_props, ANY))
        p2.assert_has_calls(exp_calls)

    @patch("q2_assembly.bowtie2.mapping._map_sample_reads")
    @patch("q2_assembly.bowtie2.mapping._gather_sample_data")
    def test_map_reads_to_contigs_single_local(self, p1, p2):
        input_reads = self.get_data_path("reads/single-end")
        input_index = self.get_data_path("indices/from_contigs")
        reads = SingleLanePerSampleSingleEndFastqDirFmt(input_reads, mode="r")
        index = Bowtie2IndexDirFmt(input_index, mode="r")

        for s in self.test_samples:
            self.test_samples[s]["rev"] = None
        p1.return_value = self.test_samples

        _map_reads_to_contigs(
            indexed_contigs=index,
            reads=reads,
            trim5=10,
            n=1,
            i="L,1,0.5",
            valid_mate_orientations="ff",
            mode="local",
            sensitivity="very-fast",
        )

        p1.assert_called_with(index, ANY, False)
        pd.testing.assert_frame_equal(
            p1.call_args[0][1], reads.manifest.view(pd.DataFrame)
        )

        exp_calls = []
        self.test_params_list[-1] += "-local"
        for s, s_props in self.test_samples.items():
            exp_calls.append(call(self.test_params_list, False, s, s_props, ANY))
        p2.assert_has_calls(exp_calls)

    def test_map_reads_to_contigs_single_parallel(self):
        input_index = self.get_data_path("indices/from_contigs")
        input_reads = get_relative_data_path(
            self.root_test_package, "formatted-reads/single-end"
        )

        index = Bowtie2IndexDirFmt(input_index, mode="r")
        index = Artifact.import_data("SampleData[SingleBowtie2Index]", index)

        reads = SingleLanePerSampleSingleEndFastqDirFmt(input_reads, mode="r")
        reads = Artifact.import_data("SampleData[SequencesWithQuality]", reads)

        with ParallelConfig():
            (out,) = self.map_reads_to_contigs.parallel(
                indexed_contigs=index,
                reads=reads,
                trim5=10,
                n=1,
                i="L,1,0.5",
                valid_mate_orientations="ff",
                mode="local",
                sensitivity="very-fast",
            )._result()

        out.validate()
        self.assertIs(out.format, BAMDirFmt)

    def test_map_reads_to_contigs_paired_parallel(self):
        input_index = self.get_data_path("indices/from_contigs")
        input_reads = get_relative_data_path(
            self.root_test_package, "formatted-reads/paired-end"
        )

        index = Bowtie2IndexDirFmt(input_index, mode="r")
        index = Artifact.import_data("SampleData[SingleBowtie2Index]", index)

        reads = SingleLanePerSamplePairedEndFastqDirFmt(input_reads, mode="r")
        reads = Artifact.import_data("SampleData[PairedEndSequencesWithQuality]", reads)

        with ParallelConfig():
            (out,) = self.map_reads_to_contigs.parallel(
                indexed_contigs=index,
                reads=reads,
                trim5=10,
                n=1,
                i="L,1,0.5",
                valid_mate_orientations="ff",
                mode="local",
                sensitivity="very-fast",
            )._result()

        out.validate()
        self.assertIs(out.format, BAMDirFmt)


if __name__ == "__main__":
    unittest.main()
