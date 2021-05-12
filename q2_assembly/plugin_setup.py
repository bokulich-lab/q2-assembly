# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Plugin, Citations)

from q2_assembly import __version__

citations = Citations.load('citations.bib', package='q2_assembly')

plugin = Plugin(
    name='assembly',
    version=__version__,
    website="https://github.com/bokulich-lab/q2-assembly",
    package='q2_assembly',
    description=(
        'QIIME 2 plugin for (meta)genome assembly and '
        'quality control thereof.'),
    short_description='QIIME 2 plugin for (meta)genome assembly.',
)
