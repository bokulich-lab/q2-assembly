# q2-assembly
![CI](https://github.com/bokulich-lab/q2-assembly/actions/workflows/ci.yml/badge.svg)

QIIME 2 plugin for (meta)genome assembly.

## Installation
Before _q2-assembly_ is available from conda, follow the installation steps described below.

```shell
curl -sLH 'Accept: application/vnd.github.v3.raw' https://api.github.com/repos/bokulich-lab/q2-assembly/contents/requirements.txt > requirements.txt
```
```shell
mamba create -yn q2-shotgun \
  -c https://packages.qiime2.org/qiime2/2022.4/tested \
  -c bioconda -c conda-forge -c default \
  --file requirements.txt

conda activate q2-shotgun
```

Additionally, a patched version of _quast_ is required, plus
q2-types-genomics and q2-assembly itself:
```shell
pip install \
  git+https://github.com/ablab/quast.git@de7e1f4891a3487f3c0df6ae27cbfba38734d686 \
  git+https://github.com/bokulich-lab/q2-types-genomics.git \
  git+https://github.com/bokulich-lab/q2-assembly.git
```

Refresh cache and check that everything worked:
```shell
qiime dev refresh-cache
qiime info
```
