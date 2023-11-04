# q2-assembly
![CI](https://github.com/bokulich-lab/q2-assembly/actions/workflows/ci-dev.yaml/badge.svg)
[![codecov](https://codecov.io/gh/bokulich-lab/q2-assembly/branch/main/graph/badge.svg?token=THMBOFUZR0)](https://codecov.io/gh/bokulich-lab/q2-assembly)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

QIIME 2 plugin for (meta)genome assembly.

## Installation
To install _q2-assembly_, follow the installation steps described below.

```shell
mamba create -yn q2-shotgun \
  -c https://packages.qiime2.org/qiime2/2022.8/tested \
  -c bioconda -c conda-forge -c default q2-assembly q2cli

conda activate q2-shotgun
```

Refresh cache and check that everything worked:
```shell
qiime dev refresh-cache
qiime info
```

## Functionality
This QIIME 2 plugin contains actions used to assemble (meta)genomes from short single/paired-end
sequencing reads. Currently, two assemblers are supported: SPAdes and MEGAHIT (for details on
the implementation and usage, please refer to the respective documentation). Below you will
find an overview of actions available in the plugin.

| Action               | Description                                                | Underlying tool                                        |
|----------------------|------------------------------------------------------------|--------------------------------------------------------|
| assemble-megahit     | Assemble contigs using MEGAHIT.                            | [MEGAHIT](https://github.com/voutcn/megahit)           |
| assemble-spades      | Assemble contigs using SPAdes.                             | [SPAdes](https://github.com/ablab/spades)              |
| evaluate-contigs     | Evaluate quality of the assembled contigs using metaQUAST. | [QUAST](https://github.com/ablab/quast)                |
| generate-reads       | Simulate NGS reads using InSilicoSeq.                      | [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq) |
| index-contigs        | Index contigs using Bowtie 2.                              | [Bowtie 2](https://github.com/BenLangmead/bowtie2)     |
| index-mags           | Index MAGs using Bowtie 2.                                 | [Bowtie 2](https://github.com/BenLangmead/bowtie2)     |
| map-reads-to-contigs | Map reads to contigs using Bowtie 2.                       | [Bowtie 2](https://github.com/BenLangmead/bowtie2)     |

## Dev environment
This repository follows the _black_ code style. To make the development slightly easier
there are a couple of pre-commit hooks included here that will ensure that your changes
follow that formatting style. Before you start working on the code, please
install the hooks by executing `make dev` in your conda environment. From then on,
they will be run automatically every time you commit any changes.
