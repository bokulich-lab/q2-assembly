# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer

setup(
    name="q2-assembly",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="BSD-3-Clause",
    packages=find_packages(),
    author="Michal Ziemski",
    author_email="ziemski.michal@gmail.com",
    description="QIIME 2 plugin for (meta)genome assembly.",
    url="https://github.com/bokulich-lab/q2-assembly",
    entry_points={"qiime2.plugins": ["q2-assembly=q2_assembly.plugin_setup:plugin"]},
    package_data={
        "q2_assembly": ["citations.bib", "assets/quast/*"],
        "q2_assembly.tests": [
            "data/*",
            "data/reads/paired-end/*",
            "data/reads/single-end/*",
            "data/html-files/*",
            "data/html-files/fake-reports/*/*",
            "data/html-files/fake-reports/*",
            "data/contigs/*",
            "data/references/*",
            "data/zip_test_data/expected/folder_0/*",
            "data/zip_test_data/expected/folder_1/subfolder_0/*",
            "data/zip_test_data/expected/folder_1/subfolder_1/*",
            "data/zip_test_data/expected/folder_1/subfolder_2/*",
            "data/zip_test_data/expected/folder_2/subfolder_0/subsubfolder_0/*",
            "data/zip_test_data/expected/folder_2/subfolder_0/subsubfolder_1/*",
            "data/zip_test_data/expected/folder_2/subfolder_0/subsubfolder_2/*",
            "data/zip_test_data/expected/folder_2/subfolder_1/subsubfolder_0/*",
            "data/zip_test_data/expected/folder_2/subfolder_1/subsubfolder_1/*",
            "data/zip_test_data/expected/folder_2/subfolder_1/subsubfolder_2/*",
            "data/zip_test_data/expected/folder_2/subfolder_2/subsubfolder_0/*",
            "data/zip_test_data/expected/folder_2/subfolder_2/subsubfolder_1/*",
            "data/zip_test_data/expected/folder_2/subfolder_2/subsubfolder_2/*",
        ],
        "q2_assembly.bowtie2.tests": [
            "data/*",
            "data/contigs/*",
            "data/indices/*/*/*",
            "data/mags/*/*",
            "data/maps/*/*",
            "data/reads/*/*",
        ],
        "q2_assembly.iss.tests": ["data/*", "data/*/*"],
    },
    zip_safe=False,
)
