import re
from functools import reduce
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer
from loguru import logger
from matplotlib.colors import ListedColormap
from typing_extensions import Annotated

BED_COLUMNS = ["chrom", "start", "end", "transcript_id"]

READS_COV = [1, 5, 10, 30, 50, 100]
LOCATION_COLS = ["CHROM", "POS", "REF", "ALT"]


def main(
    dp_file: Path,
    sample_file: Path,
    out_file_prefix: Path,
) -> None:
    sample_list = pd.read_csv(sample_file, header=None)[0].tolist()
    dp_df = pd.read_table(dp_file, header=None, names=[*LOCATION_COLS, *sample_list])
    dp_df = dp_df.set_index(LOCATION_COLS)
    dp_df.replace(".", 0, inplace=True)
    df_matrix_bool = dp_df > 0
    passed_df = df_matrix_bool.sum().reset_index()
    passed_df.columns = ["样本", "检出位点数"]
    passed_df["检出率"] = passed_df["检出位点数"] / dp_df.shape[0]
    passed_df["缺失率"] = 1 - passed_df["检出率"]
    passed_df["群体位点数"] = dp_df.shape[0]
    passed_df["缺失位点数"] = dp_df.shape[0] - passed_df["检出位点数"]
    dp_df.to_excel(f"{out_file_prefix}.目标位点深度.xlsx", index=False)
    passed_df.to_excel(
        f"{out_file_prefix}.目标位点统计.xlsx",
        index=False,
        columns=["样本", "群体位点数", "缺失位点数", "缺失率", "检出率"],
    )


if __name__ == "__main__":
    typer.run(main)
