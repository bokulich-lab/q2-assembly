# q2-assembly
![CI](https://github.com/bokulich-lab/q2-assembly/actions/workflows/ci.yaml/badge.svg)
[![codecov](https://codecov.io/gh/bokulich-lab/q2-assembly/branch/main/graph/badge.svg?token=THMBOFUZR0)](https://codecov.io/gh/bokulich-lab/q2-assembly)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

QIIME 2 plugin for (meta)genome assembly.

## Installation
_q2-assembly_ is available as part of the QIIME 2 moshpit distribution. For installation and usage instructions please consult the official [QIIME 2 documentation](https://docs.qiime2.org).

## Functionality
This QIIME 2 plugin contains actions used to assemble (meta)genomes from short single/paired-end
sequencing reads:

| Action                 | Description                                                | Underlying tool                                        |
|------------------------|------------------------------------------------------------|--------------------------------------------------------|
| assemble-megahit       | Assemble contigs using MEGAHIT.                            | [MEGAHIT](https://github.com/voutcn/megahit)           |
| assemble-spades        | Assemble contigs using SPAdes.                             | [SPAdes](https://github.com/ablab/spades)              |
| evaluate-quast         | Evaluate quality of the assembled contigs using metaQUAST. | [QUAST](https://github.com/ablab/quast)                |
| filter-contigs         | Filter contigs by length and/or metadata.                  | -                                                      |  
| generate-reads         | Simulate NGS reads using InSilicoSeq.                      | [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq) |
| index-contigs          | Index contigs using Bowtie 2.                              | [Bowtie 2](https://github.com/BenLangmead/bowtie2)     |
| index-derep-mags       | Index dereplicated MAGs using Bowtie2.                     | [Bowtie 2](https://github.com/BenLangmead/bowtie2)     |
| index-mags             | Index MAGs using Bowtie 2.                                 | [Bowtie 2](https://github.com/BenLangmead/bowtie2)     |
| map-reads              | Map reads to contigs/MAGs using Bowtie 2.                  | [Bowtie 2](https://github.com/BenLangmead/bowtie2)     |
| rename-contigs         | Rename contigs using unique IDs.                           | -                                                      |
| simulate-reads-mason   | Simulate short reads using Mason.                          | [Mason](https://www.seqan.de/apps/mason.html)          |

## Dev environment
This repository follows the _black_ code style. To make the development slightly easier
there are a couple of pre-commit hooks included here that will ensure that your changes
follow that formatting style. Before you start working on the code, please
install the hooks by executing `make dev` in your conda environment. From then on,
they will be run automatically every time you commit any changes.
