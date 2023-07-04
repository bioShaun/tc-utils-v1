from pathlib import Path
from typing import Optional
import typer

import pandas as pd
from scipy import stats


def zscoreAndFlat(
    tpm_file: Path, out_file: Path, sep="\t", db_name: Optional[str] = None
) -> None:
    tpm_df = pd.read_csv(tpm_file, sep=sep, index_col=0)
    order_cols = sorted(tpm_df.columns)
    order_tpm_df = tpm_df[order_cols]
    zscore_df = pd.DataFrame(stats.zscore(order_tpm_df.T).T)
    flat_tpm_df = order_tpm_df.melt(ignore_index=False, var_name="sample_name")
    flat_tpm_df["value_type"] = "tpm_raw"
    flat_zscore_df = zscore_df.melt(ignore_index=False, var_name="sample_name")
    flat_zscore_df["value_type"] = "tpm_zscore"
    merged_df = pd.concat([flat_tpm_df, flat_zscore_df])
    if db_name:
        merged_df["db_name"] = db_name
    merged_df.to_csv(out_file, sep=sep, float_format="%.3f")


if __name__ == "__main__":
    typer.run(zscoreAndFlat)
