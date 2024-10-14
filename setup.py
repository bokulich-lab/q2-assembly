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
            "data/reads/single-sample/paired-end/*",
            "data/reads/single-sample/single-end/*",
            "data/reads/small-single-end/*",
            "data/alignment_map/*",
            "data/html-files/*",
            "data/quast-results/*",
            "data/html-files/fake-reports/*/*",
            "data/html-files/fake-reports/*",
            "data/contigs/*",
            "data/references/*",
            "data/zip_test_data/expected/*",
            "data/zip_test_data/expected/*/*",
            "data/zip_test_data/expected/*/*/*",
            "data/zip_test_data/expected/*/*/*/*",
            "data/formatted-reads/single-end/*",
            "data/formatted-reads/paired-end/*",
            "data/dna-fasta-format/*",
            "data/genomes-dir-format1/*",
            "data/genomes-dir-format2/*",
        ],
        "q2_assembly.bowtie2.tests": [
            "data/*",
            "data/contigs/*",
            "data/empty_contigs/*",
            "data/indices/*/*/*",
            "data/indices/*/*",
            "data/mags/*/*",
            "data/mags-derep/*",
            "data/mags-merged/*",
            "data/maps/*/*",
            "data/reads/*/*",
        ],
        "q2_assembly.quast.types.tests": ["data/*"],
        "q2_assembly.iss.tests": ["data/*", "data/*/*"],
    },
    zip_safe=False,
)
