import pandas as pd

from q2_assembly.quast.report import MANDATORY_COLS_MAP, initialize_optional_cols_map
from q2_assembly.quast.types import QUASTResultsFormat


def _parse_columns(
    report_df: pd.DataFrame, contig_thresholds: list = None
) -> pd.DataFrame:
    """
    This function will rename and select the needed columns of the QUAST
    results.

    Args:
        - report_df(pd.Dataframe): Dataframe containing the QUAST results
        - contig_thresholds(list): list of contig thresholds. Defaults to None.

    Returns:
        a Pandas dataframe with the renamed columns.
    """
    report_df_newcols = report_df.copy()
    report_df_newcols.rename(columns=MANDATORY_COLS_MAP, inplace=True)
    optional_cols = []
    optional_cols_map = initialize_optional_cols_map(contig_thresholds)

    # find the optional columns
    optional_columns_present = set(optional_cols_map.keys()).intersection(
        report_df_newcols.columns
    )
    for col in optional_columns_present:
        report_df_newcols.rename(columns={col: optional_cols_map[col]}, inplace=True)
        optional_cols.append(optional_cols_map[col])

    # make sure that in the final table we have the values we
    # specify in the dicts in a certain order
    all_cols = QUASTResultsFormat.HEADER + optional_cols
    report_df_newcols = report_df_newcols[all_cols]
    report_df_newcols["id"] = report_df_newcols["id"].str.replace("_contigs", "")

    report_df_newcols = report_df_newcols.set_index("id")

    return report_df_newcols
