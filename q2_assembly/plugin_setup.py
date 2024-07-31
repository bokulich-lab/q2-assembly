# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

from q2_types.feature_data import FeatureData, Sequence
from q2_types.feature_data_mag import MAG, Contig
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.genome_data import DNASequence, GenomeData
from q2_types.per_sample_sequences import (
    Contigs,
    MAGs,
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
    SingleBowtie2Index,
)
from q2_types.per_sample_sequences._type import AlignmentMap
from q2_types.sample_data import SampleData
from qiime2.core.type import Bool, Choices, Properties, Str, TypeMap, Visualization
from qiime2.plugin import Citations, Collection, Int, List, Plugin, Range

import q2_assembly
from q2_assembly import __version__
from q2_assembly._action_params import (
    bowtie2_indexing_param_descriptions,
    bowtie2_indexing_params,
    bowtie2_mapping_param_descriptions,
    bowtie2_mapping_params,
    iss_param_descriptions,
    iss_params,
    megahit_param_descriptions,
    megahit_params,
    partition_param_descriptions,
    partition_params,
    quast_param_descriptions,
    quast_params,
    spades_param_descriptions,
    spades_params,
)
from q2_assembly.quast.types import (
    QUASTResults,
    QUASTResultsDirectoryFormat,
    QUASTResultsFormat,
)

citations = Citations.load("citations.bib", package="q2_assembly")

plugin = Plugin(
    name="assembly",
    version=__version__,
    website="https://github.com/bokulich-lab/q2-assembly",
    package="q2_assembly",
    description=(
        "QIIME 2 plugin for (meta)genome assembly and " "quality control thereof."
    ),
    short_description="QIIME 2 plugin for (meta)genome assembly.",
)

P_coassemble, T_coassembled_seqs = TypeMap(
    {
        Bool % Choices(True): FeatureData[Contig],
        Bool % Choices(False): SampleData[Contigs],
    }
)

plugin.pipelines.register_function(
    function=q2_assembly.megahit.assemble_megahit,
    inputs={"seqs": SampleData[SequencesWithQuality | PairedEndSequencesWithQuality]},
    parameters={**megahit_params, "coassemble": P_coassemble, **partition_params},
    outputs=[("contigs", T_coassembled_seqs)],
    input_descriptions={"seqs": "The paired- or single-end sequences to be assembled."},
    parameter_descriptions={
        **megahit_param_descriptions,
        **partition_param_descriptions,
    },
    output_descriptions={"contigs": "The resulting assembled contigs."},
    name="Assemble contigs using MEGAHIT.",
    description="This method uses MEGAHIT to assemble provided paired- or "
    "single-end NGS reads into contigs.",
    citations=[citations["Li2015"], citations["Li2016"]],
)

plugin.methods.register_function(
    function=q2_assembly.megahit._assemble_megahit,
    inputs={"seqs": SampleData[SequencesWithQuality | PairedEndSequencesWithQuality]},
    parameters={**megahit_params, "coassemble": P_coassemble},  # megahit_params,
    outputs=[("contigs", T_coassembled_seqs)],
    input_descriptions={"seqs": "The paired- or single-end sequences to be assembled."},
    parameter_descriptions=megahit_param_descriptions,
    output_descriptions={"contigs": "The resulting assembled contigs."},
    name="Assemble contigs using MEGAHIT.",
    description="This method uses MEGAHIT to assemble provided paired- or "
    "single-end NGS reads into contigs.",
    citations=[citations["Li2015"], citations["Li2016"]],
)

plugin.methods.register_function(
    function=q2_assembly.helpers.partition_contigs,
    inputs={"contigs": SampleData[Contigs]},
    parameters={"num_partitions": Int % Range(1, None)},
    outputs={"partitioned_contigs": Collection[SampleData[Contigs]]},
    input_descriptions={"contigs": "The contigs to partition."},
    parameter_descriptions={
        "num_partitions": "The number of partitions to split the contigs"
        " into. Defaults to partitioning into individual"
        " samples."
    },
    name="Partition contigs",
    description="Partition contigs into individual samples or the number of"
    " partitions specified.",
)

plugin.methods.register_function(
    function=q2_assembly.helpers.collate_contigs,
    inputs={"contigs": List[SampleData[Contigs]]},
    parameters={},
    outputs={"collated_contigs": SampleData[Contigs]},
    input_descriptions={"contigs": "A collection of contigs to be collated."},
    name="Collate contigs",
    description="Takes a collection of SampleData[Contigs] and collates them"
    " into a single artifact.",
)

plugin.methods.register_function(
    function=q2_assembly.helpers.rename_contigs,
    inputs={"contigs": SampleData[Contigs]},
    parameters={"uuid_type": Str % Choices(["shortuuid", "uuid3", "uuid4", "uuid5"])},
    outputs={"renamed_contigs": SampleData[Contigs]},
    input_descriptions={"contigs": "The contigs to be renamed."},
    name="Rename contigs using unique IDs.",
    description="Takes contigs for each samples in SampleData[Contigs] "
    "and renames them by changing their IDs using one of the following "
    "functions: shortuuid, uuid3, uuid4, uuid5.",
)

plugin.methods.register_function(
    function=q2_assembly.spades.assemble_spades,
    inputs={"seqs": SampleData[SequencesWithQuality | PairedEndSequencesWithQuality]},
    parameters={**spades_params, "coassemble": P_coassemble},
    outputs=[("contigs", T_coassembled_seqs)],
    input_descriptions={"seqs": "The paired- or single-end sequences to be assembled."},
    parameter_descriptions=spades_param_descriptions,
    output_descriptions={"contigs": "The resulting assembled contigs."},
    name="Assemble contigs using SPAdes.",
    description="This method uses SPAdes to assemble provided paired- or "
    "single-end NGS reads into contigs.",
    citations=[citations["Clark2021"]],
)

plugin.visualizers.register_function(
    function=q2_assembly.quast._visualize_quast,
    inputs={
        "contigs": SampleData[Contigs],
        "reads": SampleData[SequencesWithQuality | PairedEndSequencesWithQuality],
        "references": GenomeData[DNASequence],  # List[FeatureData[Sequence]],
        "mapped_reads": SampleData[AlignmentMap],
    },
    parameters=quast_params,
    input_descriptions={
        "contigs": "Assembled contigs to be analyzed.",
        "reads": "Original single- or paired-end reads.",
        "references": "Reference genomes to align the assembled contigs against.",
        "mapped_reads": "Reads-to-contigs alignment maps (alternative to 'reads')."
        "directly.",
    },
    parameter_descriptions=quast_param_descriptions,
    name="Visualize the quality of the assembled contigs after using metaQUAST.",
    description="This method visualizes the results of metaQUAST after "
    "assessing the quality of assembled metagenomes.",
    citations=[citations["Mikheenko2016"], citations["Mikheenko2018"]],
)

plugin.pipelines.register_function(
    function=q2_assembly.quast.evaluate_contigs,
    inputs={
        "contigs": SampleData[Contigs],
        "reads": SampleData[SequencesWithQuality | PairedEndSequencesWithQuality],
        "references": GenomeData[DNASequence],
        "mapped_reads": SampleData[AlignmentMap],
    },
    parameters=quast_params,
    outputs=[
        ("results_table", QUASTResults),
        ("visualization", Visualization),
        ("genomes", GenomeData[DNASequence]),
    ],
    input_descriptions={
        "contigs": "Assembled contigs to be analyzed.",
        "reads": "Original single- or paired-end reads.",
        "references": "Reference genomes to align the assembled contigs against.",
        "mapped_reads": "Reads-to-contigs alignment maps (alternative to 'reads')."
        "directly.",
    },
    parameter_descriptions=quast_param_descriptions,
    output_descriptions={
        "results_table": "QUAST result table.",
        "visualization": "Visualization of the QUAST results.",
        "genomes": "Genome sequences downloaded by QUAST.",
    },
    name="Evaluate quality of the assembled contigs using metaQUAST.",
    description="This method uses metaQUAST to assess the quality of "
    "assembled metagenomes.",
    citations=[citations["Mikheenko2016"], citations["Mikheenko2018"]],
)

plugin.pipelines.register_function(
    function=q2_assembly.indexing.index_contigs,
    inputs={"contigs": SampleData[Contigs]},
    parameters={**bowtie2_indexing_params, **partition_params},
    outputs=[("index", SampleData[SingleBowtie2Index % Properties("contigs")])],
    input_descriptions={"contigs": "Contigs to be indexed."},
    parameter_descriptions={
        **bowtie2_indexing_param_descriptions,
        **partition_param_descriptions,
    },
    output_descriptions={"index": "Bowtie2 indices generated for input sequences."},
    name="Index contigs using Bowtie2.",
    description="This method uses Bowtie2 to generate indices of " "provided contigs.",
    citations=[citations["Langmead2012"]],
)

plugin.methods.register_function(
    function=q2_assembly.indexing._index_contigs,
    inputs={"contigs": SampleData[Contigs]},
    parameters=bowtie2_indexing_params,
    outputs=[("index", SampleData[SingleBowtie2Index % Properties("contigs")])],
    input_descriptions={"contigs": "Contigs to be indexed."},
    parameter_descriptions=bowtie2_indexing_param_descriptions,
    output_descriptions={"index": "Bowtie2 indices generated for input sequences."},
    name="Index contigs using Bowtie2.",
    description="This method uses Bowtie2 to generate indices of " "provided contigs.",
    citations=[citations["Langmead2012"]],
)

I_property, O_property = TypeMap(
    {
        Properties(["mags", "contigs"]): Properties(["mags", "contigs"]),
        Properties("contigs"): Properties("contigs"),
        Properties("mags"): Properties("mags"),
    }
)
plugin.methods.register_function(
    function=q2_assembly.helpers.collate_indices,
    inputs={"indices": List[SampleData[SingleBowtie2Index % I_property]]},
    parameters={},
    outputs={"collated_indices": SampleData[SingleBowtie2Index % O_property]},
    input_descriptions={"indices": "A collection of indices to be collated."},
    name="Collate indices",
    description=(
        "Takes a collection of SampleData[Bowtie2Incex] and collates"
        " them into a single artifact."
    ),
)

plugin.methods.register_function(
    function=q2_assembly.indexing.index_mags,
    inputs={"mags": SampleData[MAGs]},
    parameters=bowtie2_indexing_params,
    outputs=[("index", SampleData[SingleBowtie2Index % Properties("mags")])],
    input_descriptions={"mags": "MAGs to be indexed."},
    parameter_descriptions=bowtie2_indexing_param_descriptions,
    output_descriptions={"index": "Bowtie2 indices generated for input sequences."},
    name="Index MAGs using Bowtie2.",
    description="This method uses Bowtie2 to generate indices of provided MAGs. One "
    "index per sample will be generated from all the MAGs belonging to that sample.",
    citations=[citations["Langmead2012"]],
)

plugin.methods.register_function(
    function=q2_assembly.indexing.index_derep_mags,
    inputs={"mags": FeatureData[MAG]},
    parameters={**bowtie2_indexing_params},
    outputs=[("index", FeatureData[SingleBowtie2Index % Properties("mags")])],
    input_descriptions={"mags": "Dereplicated MAGs to be indexed."},
    parameter_descriptions=bowtie2_indexing_param_descriptions,
    output_descriptions={"index": "Bowtie2 indices generated for input sequences."},
    name="Index dereplicated MAGs using Bowtie2.",
    description="This method uses Bowtie2 to generate indices of provided MAGs.",
    citations=[citations["Langmead2012"]],
)

plugin.methods.register_function(
    function=q2_assembly.iss.generate_reads,
    inputs={"genomes": FeatureData[Sequence]},
    parameters=iss_params,
    outputs=[
        ("reads", SampleData[PairedEndSequencesWithQuality]),
        ("template_genomes", FeatureData[Sequence]),
        ("abundances", FeatureTable[Frequency]),
    ],
    input_descriptions={
        "genomes": "Input genome(s) from which the reads will "
        "originate. If the genomes are not provided, "
        "they will be fetched from NCBI based on the "
        '"ncbi" and "n-genomes-ncbi" parameters.'
    },
    parameter_descriptions=iss_param_descriptions,
    output_descriptions={
        "reads": "Simulated paired-end reads.",
        "template_genomes": "Genome sequences from which the reads " "were generated.",
        "abundances": "Abundances of genomes from which thereads were "
        'generated. If "coverage" parameter was set, this table '
        "becomes coverage distribution per sample.",
    },
    name="Simulate NGS reads using InSilicoSeq.",
    description="This method uses InSilicoSeq to generate reads simulated "
    "from given genomes for an indicated number of samples.",
    citations=[citations["Gourle2019"]],
)

I_index, O_alignment = TypeMap(
    {
        SampleData[SingleBowtie2Index]: SampleData[AlignmentMap],
        FeatureData[SingleBowtie2Index]: FeatureData[AlignmentMap],
    }
)
plugin.pipelines.register_function(
    function=q2_assembly.mapping.map_reads,
    inputs={
        "index": I_index,
        "reads": SampleData[PairedEndSequencesWithQuality | SequencesWithQuality],
    },
    parameters={**bowtie2_mapping_params, **partition_params},
    outputs=[("alignment_map", O_alignment)],
    input_descriptions={
        "index": "Bowtie 2 indices generated for contigs/MAGs of interest.",
        "reads": "The paired- or single-end reads from which the contigs "
        "were assembled.",
    },
    parameter_descriptions={
        **bowtie2_mapping_param_descriptions,
        **partition_param_descriptions,
    },
    output_descriptions={"alignment_map": "Reads-to-contigs mapping."},
    name="Map reads to contigs using Bowtie2.",
    description="This method uses Bowtie2 to map provided reads to "
    "respective contigs.",
    citations=[citations["Langmead2012"]],
)

plugin.methods.register_function(
    function=q2_assembly.mapping._map_reads_to_contigs,
    inputs={
        "index": SampleData[SingleBowtie2Index % Properties("contigs")],
        "reads": SampleData[PairedEndSequencesWithQuality | SequencesWithQuality],
    },
    parameters=bowtie2_mapping_params,
    outputs=[("alignment_map", SampleData[AlignmentMap])],
    input_descriptions={
        "index": "Bowtie 2 indices generated for contigs of interest.",
        "reads": "The paired- or single-end reads from which the contigs "
        "were assembled.",
    },
    parameter_descriptions=bowtie2_mapping_param_descriptions,
    output_descriptions={"alignment_map": "Reads-to-contigs mapping."},
    name="Map reads to contigs using Bowtie2.",
    description="This method uses Bowtie2 to map provided reads to "
    "respective contigs.",
    citations=[citations["Langmead2012"]],
)

I_index, O_map = TypeMap(
    {
        SampleData[SingleBowtie2Index % Properties("mags")]: SampleData[AlignmentMap],
        FeatureData[SingleBowtie2Index % Properties("mags")]: FeatureData[AlignmentMap],
    }
)
plugin.methods.register_function(
    function=q2_assembly.mapping._map_reads_to_mags,
    inputs={
        "index": I_index,
        "reads": SampleData[PairedEndSequencesWithQuality | SequencesWithQuality],
    },
    parameters=bowtie2_mapping_params,
    outputs=[("alignment_map", O_map)],
    input_descriptions={
        "index": "Bowtie 2 indices generated for MAGs of interest.",
        "reads": "The paired- or single-end reads from which the contigs "
        "were assembled.",
    },
    parameter_descriptions=bowtie2_mapping_param_descriptions,
    output_descriptions={"alignment_map": "Reads-to-MAGs mapping."},
    name="Map reads to MAGs using Bowtie2.",
    description="This method uses Bowtie2 to map provided reads to "
    "the respective MAGs.",
    citations=[citations["Langmead2012"]],
)

I_maps, O_maps = TypeMap(
    {
        SampleData[AlignmentMap]: SampleData[AlignmentMap],
        FeatureData[AlignmentMap]: FeatureData[AlignmentMap],
    }
)
plugin.methods.register_function(
    function=q2_assembly.helpers.collate_alignments,
    inputs={"alignments": List[I_maps]},
    parameters={},
    outputs=[("collated_alignments", O_maps)],
    input_descriptions={"alignments": "A collection of alignments to be collated."},
    output_descriptions={
        "collated_alignments": "The alignemnts collated into one artifact."
    },
    name="Map reads to contigs helper.",
    description="Not to be called directly. Used by map_reads.",
)

plugin.register_semantic_types(QUASTResults)
plugin.register_semantic_type_to_format(
    QUASTResults, artifact_format=QUASTResultsDirectoryFormat
)
plugin.register_formats(QUASTResultsFormat, QUASTResultsDirectoryFormat)
importlib.import_module("q2_assembly.quast.types._transformer")
