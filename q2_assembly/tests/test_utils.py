# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import contextlib
import hashlib
import os
import shutil
import tempfile
import unittest
import uuid
from distutils.dir_util import copy_tree
from unittest.mock import Mock

import skbio
from bs4 import BeautifulSoup as BS
from parameterized import parameterized
from qiime2.plugin.testing import TestPluginBase

from q2_assembly._utils import (
    _construct_param,
    _generate_unique_uuid,
    _get_sample_from_path,
    _modify_links,
    _process_common_input_params,
    _remove_html_element,
    concatenate_files,
    generate_shortuuid,
    generate_uuid3,
    generate_uuid4,
    generate_uuid5,
    get_file_extension,
    modify_contig_ids,
)
from q2_assembly.quast.quast import _fix_html_reports


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

        self.namespace = uuid.NAMESPACE_OID
        self.name = "test-test"
        self.new_ids = set()

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

    def test_get_file_extension_from_path(self):

        # note that these files do not exist,
        # the file paths are used only for testing purposes
        obs1 = get_file_extension("sample_1.fastq.gz")
        exp1 = ".fastq.gz"
        obs2 = get_file_extension("sample_2.fasta.gz")
        exp2 = ".fasta.gz"
        obs3 = get_file_extension("sample_3.fastq")
        exp3 = ".fastq"
        obs4 = get_file_extension("sample_4.fa")
        exp4 = ".fa"

        self.assertEqual(obs1, exp1)
        self.assertEqual(obs2, exp2)
        self.assertEqual(obs3, exp3)
        self.assertEqual(obs4, exp4)

    def test_file_concatenation(self):
        file1 = self.get_data_path(
            "reads/small-single-end/sample1_00_L001_R1_001.fastq"
        )
        file2 = self.get_data_path(
            "reads/small-single-end/sample2_00_L001_R1_001.fastq"
        )

        # ground truth
        concat_file = self.get_data_path("reads/small-single-end/concat_samples.fastq")
        hash_original = hashlib.md5(open(concat_file, mode="rb").read()).hexdigest()

        # output of our method
        output_file = self.get_data_path(
            "reads/small-single-end/file1_file2_concat.fastq"
        )
        concatenate_files([file1, file2], output_file)
        hash_test_file = hashlib.md5(open(output_file, "rb").read()).hexdigest()

        self.assertEqual(hash_original, hash_test_file)
        os.remove(output_file)

    def test_uuid_generation(self):

        id_short = generate_shortuuid()
        id_uuid_3 = generate_uuid3(self.namespace, self.name)
        id_uuid_4 = generate_uuid4()
        id_uuid_5 = generate_uuid5(self.namespace, self.name)

        self.assertIsInstance(id_short, str)
        self.assertEqual(len(id_short), 22)

        self.assertIsInstance(id_uuid_3, uuid.UUID)
        self.assertEqual(id_uuid_3.version, 3)
        self.assertEqual(str(id_uuid_3), str(uuid.uuid3(self.namespace, self.name)))

        self.assertIsInstance(id_uuid_4, uuid.UUID)
        self.assertEqual(id_uuid_4.version, 4)

        self.assertIsInstance(id_uuid_5, uuid.UUID)
        self.assertEqual(id_uuid_5.version, 5)
        self.assertEqual(str(id_uuid_5), str(uuid.uuid5(self.namespace, self.name)))

    @parameterized.expand(["uuid3", "uuid5"])
    def test_generate_unique_uuid_uuid3_uuid5(self, uuid_type):
        mock_uuid_func = Mock(side_effect=lambda ns, cid: f"{uuid_type}-{ns}-{cid}")
        new_id = _generate_unique_uuid(
            mock_uuid_func, self.namespace, self.name, self.new_ids, uuid_type
        )
        mock_uuid_func.assert_called_with(self.namespace, self.name)
        self.assertEqual(new_id, f"{uuid_type}-{self.namespace}-{self.name}")
        self.assertIn(new_id, self.new_ids)

    @parameterized.expand(["shortuuid", "uuid4"])
    def test_generate_unique_uuid_shortuuid_uuid4(self, uuid_type):
        mock_uuid_func = Mock(side_effect=lambda: uuid_type)
        new_id = _generate_unique_uuid(
            mock_uuid_func, self.namespace, self.name, self.new_ids, uuid_type
        )
        mock_uuid_func.assert_called_once()
        self.assertEqual(new_id, uuid_type)
        self.assertIn(new_id, self.new_ids)

    def test_modify_contig_ids(self):
        contigs_path = self.get_data_path("contigs")
        new_ids_sample1 = []
        new_ids_sample2 = []
        with tempfile.TemporaryDirectory() as tmp:
            for sample in ["sample1", "sample2"]:
                original_sample_path = os.path.join(
                    contigs_path, f"{sample}_contigs.fa"
                )
                new_sample_path = os.path.join(tmp, f"{sample}_contigs.fa")
                shutil.copy(original_sample_path, tmp)
                modify_contig_ids(tmp, sample, "uuid4")

                for contig in skbio.io.read(new_sample_path, format="fasta"):
                    if sample == "sample1":
                        new_ids_sample1.append(contig.metadata["id"])
                    else:
                        new_ids_sample2.append(contig.metadata["id"])

            new_ids_sample1_set = set(new_ids_sample1)
            new_ids_sample2_set = set(new_ids_sample2)
            self.assertEqual(len(new_ids_sample1), len(new_ids_sample1_set))
            self.assertEqual(len(new_ids_sample2), len(new_ids_sample2_set))
            self.assertEqual(
                len(new_ids_sample1_set.intersection(new_ids_sample2_set)), 0
            )

    def test_generate_unique_uuid_with_existing_ids(self):
        self.new_ids.add("uuid4")
        mock_uuid_func = Mock(side_effect=["uuid4", "uuid4", "uuid4", "unique-uuid4"])
        new_id = _generate_unique_uuid(
            mock_uuid_func, self.namespace, self.name, self.new_ids, "uuid4"
        )
        self.assertEqual(mock_uuid_func.call_count, 4)
        self.assertEqual(new_id, "unique-uuid4")
        self.assertIn(new_id, self.new_ids)

    def test_get_sample_from_path(self):
        obs = _get_sample_from_path("test/path/sample_1_contigs.fa")
        exp = "sample_1"
        self.assertEqual(obs, exp)


if __name__ == "__main__":
    unittest.main()
