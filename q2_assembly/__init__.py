# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from .bowtie2 import indexing, mapping
from .filter import filter
from .helpers import helpers
from .iss import iss
from .megahit import megahit
from .quast import quast
from .spades import spades

try:
    from ._version import __version__
except ModuleNotFoundError:
    __version__ = '0.0.0+notfound'

__all__ = [
    "indexing",
    "mapping",
    "iss",
    "megahit",
    "quast",
    "spades",
    "helpers",
    "filter",
]
