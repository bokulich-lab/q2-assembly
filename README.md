# q2-assembly
![CI](https://github.com/bokulich-lab/q2-assembly/actions/workflows/ci.yml/badge.svg)
[![codecov](https://codecov.io/gh/bokulich-lab/q2-assembly/branch/main/graph/badge.svg?token=THMBOFUZR0)](https://codecov.io/gh/bokulich-lab/q2-assembly)

QIIME 2 plugin for (meta)genome assembly.

## Installation
To install _q2-assembly_, follow the installation steps described below.

```shell
mamba create -yn q2-shotgun \
  -c https://packages.qiime2.org/qiime2/2022.4/tested \
  -c bioconda -c conda-forge -c default q2-assembly q2cli

conda activate q2-shotgun
```

Refresh cache and check that everything worked:
```shell
qiime dev refresh-cache
qiime info
```
