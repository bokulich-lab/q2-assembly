# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from .._utils import _construct_param, _process_common_input_params


def fake_processing_func(key, val):
    return [_construct_param(key), str(val)]


class TestUtils(unittest.TestCase):
    def test_construct_param_simple(self):
        obs = _construct_param('test')
        exp = '--test'
        self.assertEqual(obs, exp)

    def test_construct_param_complex(self):
        obs = _construct_param('test_param')
        exp = '--test-param'
        self.assertEqual(obs, exp)

    def test_process_common_inputs_bools(self):
        kwargs = {'arg1': False, 'arg2': True}
        obs = _process_common_input_params(fake_processing_func, **kwargs)
        exp = ['--arg2']
        self.assertListEqual(obs, exp)

    def test_process_common_inputs_nones(self):
        kwargs = {'arg1': 'some-value', 'arg2': None}
        obs = _process_common_input_params(fake_processing_func, **kwargs)
        exp = ['--arg1', 'some-value']
        self.assertListEqual(obs, exp)

    def test_process_common_inputs_with_values(self):
        kwargs = {'arg1': 'value1', 'arg2': 'value2'}
        obs = _process_common_input_params(fake_processing_func, **kwargs)
        exp = ['--arg1', 'value1', '--arg2', 'value2']
        self.assertListEqual(obs, exp)

    def test_process_common_inputs_mix(self):
        kwargs = {
            'arg1': None, 'arg2': 'some-value', 'arg3': False, 'arg4': True
        }
        obs = _process_common_input_params(fake_processing_func, **kwargs)
        exp = ['--arg2', 'some-value', '--arg4']
        self.assertListEqual(obs, exp)


if __name__ == '__main__':
    unittest.main()
