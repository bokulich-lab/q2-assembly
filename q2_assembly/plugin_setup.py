# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.per_sample_sequences import (SequencesWithQuality,
                                           PairedEndSequencesWithQuality)
from q2_types.sample_data import SampleData
from q2_types_genomics.per_sample_data import Contigs
from qiime2.core.type import Str, Int, List, Range, Bool, Float
from qiime2.plugin import (Plugin, Citations)

import q2_assembly
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

plugin.methods.register_function(
    function=q2_assembly.megahit.assemble_megahit,
    inputs={
        'seqs': SampleData[SequencesWithQuality |
                           PairedEndSequencesWithQuality]
    },
    parameters={
        'presets': Str,
        'min_count': Int % Range(1, None),
        'k_list': List[Int % Range(15, 255, inclusive_end=True)],
        'k_min': Int % Range(15, 255, inclusive_end=True),
        'k_max': Int % Range(15, 255, inclusive_end=True),
        'k_step': Int % Range(2, 28, inclusive_end=True),
        'no_mercy': Bool,
        'bubble_level': Int % Range(0, 2, inclusive_end=True),
        'prune_level': Int % Range(0, 3, inclusive_end=True),
        'prune_depth': Int % Range(1, None),
        'disconnect_ratio': Float % Range(0, 1, inclusive_end=True),
        'low_local_ratio': Float % Range(0, 1, inclusive_end=True),
        'max_tip_len': Int % Range(1, None),
        'cleaning_rounds': Int % Range(1, None),
        'no_local': Bool,
        'kmin_1pass': Bool,
        'memory': Float % Range(0, None),
        'mem_flag': Int % Range(0, None),
        'num_cpu_threads': Int % Range(1, None),
        'no_hw_accel': Bool,
        'min_contig_len': Int
    },
    outputs=[('contigs', SampleData[Contigs])],
    input_descriptions={'seqs': 'The paired- or single-end sequences '
                                'to be assembled.'},
    parameter_descriptions={
        'presets': 'Override a group of parameters. Possible values: '
                   '"meta-sensitive", "meta-large".',
        'min_count': 'Minimum multiplicity for filtering (k_min+1)-mers. '
                     'Default: 2',
        'k_list': 'List of kmer size - all must be odd with an increment '
                  '<= 28. Default: [21,29,39,59,79,99,119,141]',
        'k_min': 'Minimum kmer size (<= 255), must be odd number. '
                 'Default: 21. Overrides k_list.',
        'k_max': 'Maximum kmer size (<= 255), must be odd number. '
                 'Default: 141. Overrides k_list.',
        'k_step': 'Increment of kmer size of each iteration (<= 28), '
                  'must be even number. Default: 12. Overrides k_list.',
        'no_mercy': 'Do not add mercy kmers.',
        'bubble_level': 'Intensity of bubble merging, 0 to disable. '
                        'Default: 2.',
        'prune_level': 'Strength of low depth pruning. Default: 2.',
        'prune_depth': 'Remove unitigs with avg kmer depth less than '
                       'this value. Default: 2.',
        'disconnect_ratio': 'Disconnect unitigs if its depth is less '
                            'than this ratio times the total depth of '
                            'itself and its siblings. Default: 0.1.',
        'low_local_ratio': 'Remove unitigs if its depth is less than '
                           'this ratio times the average depth of '
                           'the neighborhoods. Default: 0.2.',
        'max_tip_len': 'Remove tips less than this value.',
        'cleaning_rounds': 'Number of rounds for graph cleanning. Default: 5.',
        'no_local': 'Disable local assembly.',
        'kmin_1pass': 'Use 1pass mode to build SdBG of k_min.',
        'memory': 'Max memory in byte to be used in SdBG construction '
                  '(if set between 0-1, fraction of the machine\'s total '
                  'memory). Default: 0.9.',
        'mem_flag': 'SdBG builder memory mode. 0: minimum; 1: moderate; '
                    'others: use all memory specified by \'-m/--memory\'. '
                    'Default: 1.',
        'num_cpu_threads': 'Number of CPU threads. '
                           'Default: # of logical processors.',
        'no_hw_accel': 'Run MEGAHIT without BMI2 and POPCNT '
                       'hardware instructions.',
        'min_contig_len': 'Minimum length of contigs to output. Default: 200.'
    },
    output_descriptions={'contigs': 'The resulting assembled contigs.'},
    name='Assemble contigs using MEGAHIT.',
    description='This method uses MEGAHIT to assemble provided paired- or '
                'single-end NGS reads into contigs.',
    citations=[
        citations['Li2015'],
        citations['Li2016']]
)
