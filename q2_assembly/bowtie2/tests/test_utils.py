# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import filecmp
import unittest

from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt
from qiime2.plugin.testing import TestPluginBase

from q2_assembly.bowtie2.utils import (
    _construct_double_list_param_value,
    _construct_function_param_value,
    _get_subdir_from_path,
    _is_flat_dir,
    _merge_mags,
    _process_bowtie2_arg,
    _process_bowtie2build_arg,
)


class TestBowtie2Utils(TestPluginBase):
    package = "q2_assembly.bowtie2.tests"

    def test_get_subdir_from_path_contig(self):
        obs = _get_subdir_from_path("/path/to/dir/sample1_contigs.fa")
        exp = "sample1"
        self.assertEqual(obs, exp)

    def test_get_subdir_from_path_contig_underscores(self):
        self.assertEqual(
            _get_subdir_from_path("/path/to/dir/sample_1_contigs.fa"), "sample_1"
        )
        self.assertEqual(
            _get_subdir_from_path("/path/to/dir/sample1_contigs.fa"), "sample1"
        )
        self.assertEqual(
            _get_subdir_from_path("/path/to/dir/s_am-p_le_1_contigs.fa"), "s_am-p_le_1"
        )
        self.assertEqual(
            _get_subdir_from_path("/path/to/dir/_contigs.fa_contigs.fa"), "_contigs.fa"
        )

    def test_get_subdir_from_path_mag(self):
        obs = _get_subdir_from_path("/path/to/dir/sample1/mag1.fa", "mags")
        exp = "sample1"
        self.assertEqual(obs, exp)

    def test_get_subdir_from_path_mag_derep(self):
        obs = _get_subdir_from_path("/path/to/dir/mag1.fa", "mags-derep")
        exp = ""
        self.assertEqual(obs, exp)

    def test_get_subdir_from_path_unknown_type(self):
        with self.assertRaisesRegex(NotImplementedError, r'"unicorn" is not supported'):
            _ = _get_subdir_from_path("/path/to/dir/mag1.fa", "unicorn")

    def test_process_bowtie2build_arg_simple1(self):
        obs = _process_bowtie2build_arg("not_k_list", 123)
        exp = ["--not-k-list", "123"]
        self.assertListEqual(obs, exp)

    def test_process_bowtie2build_arg_simple2(self):
        with self.assertRaisesRegex(
            Exception, r'.*type "\<class \'list\'\>" is not supported\.'
        ):
            _process_bowtie2build_arg("k_list", [1, 2, 3])

    def test_process_bowtie2build_arg_bool(self):
        obs = _process_bowtie2build_arg("k_bool", True)
        exp = ["--k-bool"]
        self.assertListEqual(obs, exp)

    def test_process_bowtie2_arg_N(self):
        obs = _process_bowtie2_arg("n", 1)
        exp = ["-N", "1"]
        self.assertListEqual(obs, exp)

    def test_process_bowtie2_arg_k(self):
        obs = _process_bowtie2_arg("k", 10)
        exp = ["-k", "10"]
        self.assertListEqual(obs, exp)

    def test_process_bowtie2_arg_L(self):
        obs = _process_bowtie2_arg("len", 20)
        exp = ["-L", "20"]
        self.assertListEqual(obs, exp)

    def test_process_bowtie2_arg_a(self):
        obs = _process_bowtie2_arg("a", True)
        exp = ["-a"]
        self.assertListEqual(obs, exp)

    def test_process_bowtie2_arg_bool(self):
        obs = _process_bowtie2_arg("very_fast_local", True)
        exp = ["--very-fast-local"]
        self.assertListEqual(obs, exp)

    def test_process_bowtie2_arg_rdg(self):
        obs = _process_bowtie2_arg("rdg", "5, 3")
        exp = ["--rdg", "5,3"]
        self.assertListEqual(obs, exp)

    def test_process_bowtie2_arg_i(self):
        obs = _process_bowtie2_arg("i", "L,1,-0.5")
        exp = ["-i", "L,1,-0.5"]
        self.assertListEqual(obs, exp)

    def test_process_bowtie2_arg_n_ceil(self):
        obs = _process_bowtie2_arg("n_ceil", "L,1,-0.5")
        exp = ["--n-ceil", "L,1,-0.5"]
        self.assertListEqual(obs, exp)

    def test_process_bowtie2_arg_valid_mates(self):
        obs = _process_bowtie2_arg("valid_mate_orientations", "ff")
        exp = ["--ff"]
        self.assertListEqual(obs, exp)

    def test_construct_double_list_param_value(self):
        obs = _construct_double_list_param_value("some_param", "1,2")
        exp = "1,2"
        self.assertEqual(obs, exp)

    def test_construct_double_list_param_value_whitespaces(self):
        obs = _construct_double_list_param_value("some_param", "1, 2 ")
        exp = "1,2"
        self.assertEqual(obs, exp)

    def test_construct_double_list_param_value_too_many(self):
        with self.assertRaisesRegex(
            Exception, 'Invalid number of elements for "some_param"'
        ):
            _construct_double_list_param_value("some_param", "1,2,3")

    def test_construct_double_list_param_value_floats(self):
        with self.assertRaisesRegex(
            Exception,
            'Both values of "some_param" parameter should be '
            "integers. Provided values were: 1,0.5.",
        ):
            _construct_double_list_param_value("some_param", "1,0.5")

    def test_construct_function_param_value(self):
        obs = _construct_function_param_value("some_param", "L,0.1,-1")
        exp = "L,0.1,-1"
        self.assertEqual(obs, exp)

    def test_construct_function_param_value_whitespaces(self):
        obs = _construct_function_param_value("some_param", " L,0.1 ,-1 ")
        exp = "L,0.1,-1"
        self.assertEqual(obs, exp)

    def test_construct_function_param_value_too_many(self):
        with self.assertRaisesRegex(
            Exception, 'Invalid number .* definition of "some_param".'
        ):
            _construct_function_param_value("some_param", "A,0.1,-1,1")

    def test_construct_function_param_value_invalid_function(self):
        with self.assertRaisesRegex(
            Exception,
            'Invalid function type in "some_param": A was given '
            'but only "CLSG" are allowed.',
        ):
            _construct_function_param_value("some_param", "A,0.1,-1")

    def test_merge_mags(self):
        mags = MultiMAGSequencesDirFmt(self.get_data_path("mags"), "r")

        obs_fps = _merge_mags(mags, self.temp_dir.name)

        self.assertListEqual(
            obs_fps,
            [
                f"{self.temp_dir.name}/sample1/merged.fasta",
                f"{self.temp_dir.name}/sample2/merged.fasta",
            ],
        )
        self.assertTrue(
            filecmp.cmp(
                obs_fps[0], self.get_data_path("mags-merged/sample1.fa"), shallow=False
            )
        )
        self.assertTrue(
            filecmp.cmp(
                obs_fps[1], self.get_data_path("mags-merged/sample2.fa"), shallow=False
            )
        )

    def test_merge_mags_derep(self):
        mags = MAGSequencesDirFmt(self.get_data_path("mags-derep"), "r")

        obs_fps = _merge_mags(mags, self.temp_dir.name)

        self.assertEqual(obs_fps[0], f"{self.temp_dir.name}/merged.fasta")
        self.assertTrue(
            filecmp.cmp(
                obs_fps[0], self.get_data_path("mags-derep-merged.fasta"), shallow=False
            )
        )

    def test_is_flat_dir_false(self):
        self.assertFalse(_is_flat_dir(self.get_data_path("indices/from_mags")))

    def test_is_flat_dir_true(self):
        self.assertTrue(_is_flat_dir(self.get_data_path("indices/from_mags_derep")))


if __name__ == "__main__":
    unittest.main()
