# q2-assembly
![CI](https://github.com/bokulich-lab/q2-assembly/actions/workflows/ci.yml/badge.svg)

QIIME 2 plugin for (meta)genome assembly.

## Installation
Before _q2-assembly_ is available from conda, follow the installation steps described below.

```shell
conda create -yn q2-shotgun-env python=3.8
conda activate q2-shotgun-env
```

Some actions require updated qiime2, which you should install from here before it lands in Q2 2021.8:
```shell
pip install git+git://github.com/misialq/qiime2.git@7c8c1042b0f5dd1b6ce5e68f89541f915880f89f
```

Install remaining dependencies:
```shell
conda install -y \
  -c https://packages.qiime2.org/qiime2/2021.8/staged -c bioconda \
  -c conda-forge -c default \
  q2cli q2-types q2templates megahit beautifulsoup4 spades
pip install \
  git+git://github.com/ablab/quast.git@bc4af762a7f53176d66bd5e6c5b7d28376d28e11 \
  git+git://github.com/bokulich-lab/q2-types-genomics.git
```

Finally, install q2-assembly:
```shell
pip install git+git://github.com/bokulich-lab/q2-assembly.git
```

Refresh cache and check that everything worked:
```shell
qiime dev refresh-cache
qiime info
```
