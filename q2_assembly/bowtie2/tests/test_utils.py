# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from qiime2.plugin.testing import TestPluginBase

from q2_assembly.bowtie2.utils import (_process_bowtie2build_arg,
                                       _get_subdir_from_path)


class TestBowtie2Utils(TestPluginBase):
    package = 'q2_assembly.bowtie2.tests'

    def setUp(self):
        super().setUp()
        self.test_params_list = ['--large-index', '--bmax', '11']

    def test_get_subdir_from_path_contig(self):
        obs = _get_subdir_from_path('/path/to/dir/sample1_contigs.fa')
        exp = 'sample1'
        self.assertEqual(obs, exp)

    def test_get_subdir_from_path_mag(self):
        obs = _get_subdir_from_path('/path/to/dir/sample1/mag1.fa', 'mags')
        exp = 'sample1/mag1'
        self.assertEqual(obs, exp)

    def test_get_subdir_from_path_unknown_type(self):
        with self.assertRaisesRegex(
                NotImplementedError, r'"unicorn" is not supported'):
            _ = _get_subdir_from_path('/path/to/dir/mag1.fa', 'unicorn')

    def test_process_bowtie2_arg_simple1(self):
        obs = _process_bowtie2build_arg('not_k_list', 123)
        exp = ['--not-k-list', '123']
        self.assertListEqual(obs, exp)

    def test_process_bowtie2_arg_simple2(self):
        with self.assertRaisesRegex(
                Exception, r'.*type "\<class \'list\'\>" is not supported\.'):
            _process_bowtie2build_arg('k_list', [1, 2, 3])

    def test_process_bowtie2_arg_bool(self):
        obs = _process_bowtie2build_arg('k_bool', True)
        exp = ['--k-bool']
        self.assertListEqual(obs, exp)


if __name__ == '__main__':
    unittest.main()
