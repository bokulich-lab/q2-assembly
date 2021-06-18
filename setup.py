# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer

setup(
    name='q2-assembly',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='BSD-3-Clause',
    packages=find_packages(),
    author="Michal Ziemski",
    author_email="ziemski.michal@gmail.com",
    description="QIIME 2 plugin for (meta)genome assembly.",
    url="https://github.com/bokulich-lab/q2-assembly",
    entry_points={
        'qiime2.plugins':
        ['q2-assembly=q2_assembly.plugin_setup:plugin']
    },
    package_data={
        'q2_assembly': [
            'citations.bib'
        ],
        'q2_assembly.tests': [
            'data/*', 'data/reads/paired-end/*', 'data/reads/single-end/*'
        ]
    },
    zip_safe=False,
)
