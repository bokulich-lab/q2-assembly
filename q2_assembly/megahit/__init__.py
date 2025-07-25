# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from .megahit import _assemble_megahit, assemble_megahit

__all__ = ["assemble_megahit", "_assemble_megahit"]
