# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import csv

from qiime2.core.exceptions import ValidationError
from qiime2.plugin import model


class QUASTResultsFormat(model.TextFileFormat):
    HEADER = [
        "id",
        "total_length",
        "no_contigs_0",
        "no_contigs",
        "longest_contig",
        "n50",
        "l50",
        "n90",
        "l90",
    ]

    def _validate(self, n_records=None):
        with self.open() as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
            n_columns = len(header)

            if not set(header).issuperset(set(self.HEADER)):
                raise ValidationError(
                    f"Invalid header: {header}, do not contain "
                    f"all headers in: {self.HEADER}"
                )

            for i, row in enumerate(reader, start=2):
                if len(row) != n_columns:
                    raise ValidationError(
                        f"Line {i} has {len(row)} columns, " f"expected {n_columns}"
                    )

                if n_records is not None and i - 1 >= n_records:
                    break

    def _validate_(self, level):
        record_count_map = {"min": 10, "max": None}
        self._validate(record_count_map[level])


QUASTResultsDirectoryFormat = model.SingleFileDirectoryFormat(
    "QUASTResultsDirectoryFormat", "quast_results.tsv", QUASTResultsFormat
)
