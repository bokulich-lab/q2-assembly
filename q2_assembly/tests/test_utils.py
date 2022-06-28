# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import contextlib
import os
import shutil
import tempfile
import unittest
from distutils.dir_util import copy_tree

from bs4 import BeautifulSoup as BS
from qiime2.plugin.testing import TestPluginBase

from q2_assembly.quast.quast import _fix_html_reports

from .._utils import (
    _construct_param,
    _modify_links,
    _process_common_input_params,
    _remove_html_element,
)


def fake_processing_func(key, val):
    if not val:
        return
    elif isinstance(val, bool):
        return [_construct_param(key)]
    else:
        return [_construct_param(key), str(val)]


class TestUtils(TestPluginBase):
    package = "q2_assembly.tests"

    def setUp(self):
        super().setUp()
        with contextlib.ExitStack() as stack:
            self._tmp = stack.enter_context(tempfile.TemporaryDirectory())
            self.addCleanup(stack.pop_all().close)

    def assertHTMLEqual(self, obs_fp, exp_fp):
        with open(obs_fp, "r") as fobs, open(exp_fp, "r") as fexp:
            soup_obs = BS(fobs.read(), "html.parser")
            soup_exp = BS(fexp.read(), "html.parser")
            self.assertEqual(str(soup_obs), str(soup_exp))

    def clone_input_html(self, filepath):
        input_fp = self.get_data_path(filepath)
        shutil.copy(input_fp, self._tmp)
        output_fp = os.path.join(self._tmp, os.path.basename(input_fp))
        return input_fp, output_fp

    def test_construct_param_simple(self):
        obs = _construct_param("test")
        exp = "--test"
        self.assertEqual(obs, exp)

    def test_construct_param_complex(self):
        obs = _construct_param("test_param")
        exp = "--test-param"
        self.assertEqual(obs, exp)

    def test_process_common_inputs_bools(self):
        kwargs = {"arg1": False, "arg2": True}
        obs = _process_common_input_params(fake_processing_func, kwargs)
        exp = ["--arg2"]
        self.assertListEqual(obs, exp)

    def test_process_common_inputs_nones(self):
        kwargs = {"arg1": "some-value", "arg2": None}
        obs = _process_common_input_params(fake_processing_func, kwargs)
        exp = ["--arg1", "some-value"]
        self.assertListEqual(obs, exp)

    def test_process_common_inputs_with_values(self):
        kwargs = {"arg1": "value1", "arg2": "value2"}
        obs = _process_common_input_params(fake_processing_func, kwargs)
        exp = ["--arg1", "value1", "--arg2", "value2"]
        self.assertListEqual(obs, exp)

    def test_process_common_inputs_mix(self):
        kwargs = {"arg1": None, "arg2": "some-value", "arg3": False, "arg4": True}
        obs = _process_common_input_params(fake_processing_func, kwargs)
        exp = ["--arg2", "some-value", "--arg4"]
        self.assertListEqual(obs, exp)

    def test_remove_html_element_existing(self):
        _, fp_out = self.clone_input_html("html-files/sample-report.html")
        _remove_html_element(fp_out, "div", elem_id="div1")
        exp_fp = self.get_data_path("html-files/sample-report-div-removed.html")
        self.assertHTMLEqual(fp_out, exp_fp)

    def test_remove_html_element_nonexisting(self):
        fp_in, fp_out = self.clone_input_html("html-files/sample-report.html")
        _remove_html_element(fp_out, "div", elem_id="div8")
        self.assertHTMLEqual(fp_out, fp_in)

    def test_modify_links(self):
        fp_in, fp_out = self.clone_input_html("html-files/sample-report.html")
        _modify_links(fp_out)
        exp_fp = self.get_data_path("html-files/sample-report-links-updated.html")
        self.assertHTMLEqual(fp_out, exp_fp)

    def test_fix_html_reports(self):
        input_dir = self.get_data_path("html-files/fake-reports")
        copy_tree(input_dir, self._tmp)
        _fix_html_reports(self._tmp)
        self.assertHTMLEqual(
            f"{self._tmp}/report.html", f"{self._tmp}/expected/report.html"
        )
        self.assertHTMLEqual(
            f"{self._tmp}/icarus_viewers/contig_size_viewer.html",
            f"{self._tmp}/expected/contig_size_viewer.html",
        )


if __name__ == "__main__":
    unittest.main()
