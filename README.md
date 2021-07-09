# q2-assembly
![CI](https://github.com/bokulich-lab/q2-assembly/actions/workflows/ci.yml/badge.svg)

QIIME 2 plugin for (meta)genome assembly.

## Installation
Before _q2-assembly_ is available from conda, follow the installation steps described below.

```shell
conda create -yn q2-shotgun \
  -c https://packages.qiime2.org/qiime2/2021.8/staged -c bioconda \
  -c conda-forge -c default \
  q2cli q2-types q2templates "megahit==1.2.9" beautifulsoup4 "spades==3.15.2" \
  "bowtie2==2.4.4" insilicoseq "biopython<=1.78" samtools
conda activate q2-shotgun
```

Some actions require updated qiime2 which you should install from the fork below before 
it lands in Q2 2021.8. Additionally, a patched version of _quast_ is required, plus
q2-types-genomics and q2-assembly itself:
```shell
pip install \
  git+https://github.com/misialq/qiime2.git@7c8c1042b0f5dd1b6ce5e68f89541f915880f89f \
  git+https://github.com/ablab/quast.git@bc4af762a7f53176d66bd5e6c5b7d28376d28e11 \
  git+https://github.com/bokulich-lab/q2-types-genomics.git \
  git+https://github.com/bokulich-lab/q2-assembly.git
```

Refresh cache and check that everything worked:
```shell
qiime dev refresh-cache
qiime info
```
