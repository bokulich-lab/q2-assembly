# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from qiime2.plugin.testing import TestPluginBase


class TestBowtie2Mapping(TestPluginBase):
    package = 'q2_assembly.bowtie2.tests'

    def setUp(self):
        super().setUp()
        self.test_params_list = ['--large-index', '--bmax', '11']


if __name__ == '__main__':
    unittest.main()
