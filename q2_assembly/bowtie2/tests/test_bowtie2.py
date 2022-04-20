# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from subprocess import CalledProcessError
from unittest.mock import patch, ANY, call

from q2_types_genomics.per_sample_data import (ContigSequencesDirFmt,
                                               MultiMAGSequencesDirFmt)
from qiime2.plugin.testing import TestPluginBase

from q2_assembly.bowtie2.bowtie2 import (_get_subdir_from_path, _index_seqs,
                                         index_contigs, index_mags,
                                         _process_bowtie2_arg)


class TestBowtie2(TestPluginBase):
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
        obs = _process_bowtie2_arg('not_k_list', 123)
        exp = ['--not-k-list', '123']
        self.assertListEqual(obs, exp)

    def test_process_bowtie2_arg_simple2(self):
        with self.assertRaisesRegex(
                Exception, r'.*type "\<class \'list\'\>" is not supported\.'):
            _process_bowtie2_arg('k_list', [1, 2, 3])

    def test_process_bowtie2_arg_bool(self):
        obs = _process_bowtie2_arg('k_bool', True)
        exp = ['--k-bool']
        self.assertListEqual(obs, exp)

    @patch('subprocess.run')
    @patch('os.makedirs')
    def test_index_seqs_contigs(self, p1, p2):
        _index_seqs(
            fasta_fps=['/here/samp1_contigs.fa', '/here/samp2_contigs.fa'],
            result_fp='/there/',
            common_args=self.test_params_list,
            input_type='contigs'
        )

        p1.assert_has_calls([call('/there/samp1'), call('/there/samp2')])
        p2.assert_has_calls([
            call(['bowtie2-build', '--large-index', '--bmax', '11',
                  '/here/samp1_contigs.fa', '/there/samp1/index'], check=True),
            call(['bowtie2-build', '--large-index', '--bmax', '11',
                  '/here/samp2_contigs.fa', '/there/samp2/index'], check=True)
        ])

    @patch('subprocess.run')
    @patch('os.makedirs')
    def test_index_seqs_mags(self, p1, p2):
        _index_seqs(
            fasta_fps=['/here/smp1/mag1.fa', '/here/smp1/mag2.fa'],
            result_fp='/there/',
            common_args=self.test_params_list,
            input_type='mags'
        )

        p1.assert_has_calls([
            call('/there/smp1/mag1'), call('/there/smp1/mag2')])
        p2.assert_has_calls([
            call(['bowtie2-build', '--large-index', '--bmax', '11',
                  '/here/smp1/mag1.fa', '/there/smp1/mag1/index'], check=True),
            call(['bowtie2-build', '--large-index', '--bmax', '11',
                  '/here/smp1/mag2.fa', '/there/smp1/mag2/index'], check=True)
        ])

    @patch('subprocess.run',
           side_effect=CalledProcessError(returncode=123, cmd="some cmd"))
    @patch('os.makedirs')
    def test_index_seqs_with_error(self, p1, p2):
        with self.assertRaisesRegex(
                Exception, 'An error.*while running Bowtie2.*code 123'):
            _index_seqs(
                fasta_fps=['/here/samp1/mag1.fa'],
                result_fp='/there/',
                common_args=self.test_params_list,
                input_type='mags'
            )

    @patch('q2_assembly.bowtie2._index_seqs')
    def test_index_contigs(self, p):
        input_contigs = ContigSequencesDirFmt(
            self.get_data_path('contigs'), 'r')
        index_contigs(input_contigs, large_index=True, bmax=11)

        exp_contigs = [
            f'{str(input_contigs)}/sample{x+1}_contigs.fa' for x in range(2)]
        p.assert_called_with(
            exp_contigs, ANY, self.test_params_list, 'contigs')

    @patch('q2_assembly.bowtie2._index_seqs')
    def test_index_mags(self, p):
        input_mags = MultiMAGSequencesDirFmt(
            self.get_data_path('mags'), 'r')
        index_mags(input_mags, large_index=True, bmax=11)

        exp_mags = [
            f'{str(input_mags)}/sample{x+1}/mag{y+1}.fa'
            for x in range(2) for y in range(2)]
        p.assert_called_with(
            exp_mags, ANY, self.test_params_list, 'mags')


if __name__ == '__main__':
    unittest.main()
