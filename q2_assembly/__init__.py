# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._version import get_versions
from .bowtie2 import indexing, mapping
from .helpers import helpers
from .iss import iss
from .mason import mason
from .megahit import megahit
from .quast import quast
from .spades import spades

__version__ = get_versions()["version"]
del get_versions

__all__ = [
    "indexing",
    "mapping",
    "iss",
    "mason",
    "megahit",
    "quast",
    "spades",
    "helpers",
]
