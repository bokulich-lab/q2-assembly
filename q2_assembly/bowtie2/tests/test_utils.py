# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from qiime2.plugin.testing import TestPluginBase

from q2_assembly.bowtie2.utils import (
    _construct_double_list_param_value,
    _construct_function_param_value,
    _get_subdir_from_path,
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
        obs = _get_subdir_from_path("/path/to/dir/sample_1_contigs.fa")
        exp = "sample_1"
        self.assertEqual(obs, exp)

    def test_get_subdir_from_path_mag(self):
        obs = _get_subdir_from_path("/path/to/dir/sample1/mag1.fa", "mags")
        exp = "sample1/mag1"
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


if __name__ == "__main__":
    unittest.main()
