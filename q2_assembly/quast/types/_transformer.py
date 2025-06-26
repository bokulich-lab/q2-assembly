# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
from qiime2 import Metadata

from ...plugin_setup import plugin
from . import QUASTResultsFormat


def _read_dataframe(fh: str, header=0):
    df = pd.read_csv(fh, sep="\t", header=header, index_col=0, dtype="str")
    df.index.name = "id"
    return df


@plugin.register_transformer
def _1(ff: QUASTResultsFormat) -> pd.DataFrame:
    with ff.open() as fh:
        df = _read_dataframe(fh)
        return df


@plugin.register_transformer
def _2(data: pd.DataFrame) -> QUASTResultsFormat:
    ff = QUASTResultsFormat()
    with ff.open() as fh:
        data.to_csv(fh, sep="\t", index=False, header=True)
    return ff


@plugin.register_transformer
def _3(ff: QUASTResultsFormat) -> Metadata:
    with ff.open() as fh:
        df = _read_dataframe(fh)
        return Metadata(df)
