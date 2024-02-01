import typer
from pathlib import Path
import pandas as pd
from functools import reduce
from typing import List, Optional, Tuple
from loguru import logger
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np

BED_COLUMNS = ["chrom", "start", "end", "transcript_id"]

READS_COV = [1, 5, 10, 30, 50, 100]


def merge_chr(df: pd.DataFrame, split_bed: Path) -> pd.DataFrame:
    split_bed_df = pd.read_csv(
        split_bed,
        header=None,
        names=["new_chrom", "offset", "offset_end", "chrom"],
        sep="\t",
    )
    merged_df = df.merge(split_bed_df)
    merged_df["new_start"] = merged_df["start"] + merged_df["offset"]
    merged_df["new_end"] = merged_df["end"] + merged_df["offset"]
    merged_df.drop(
        ["chrom", "offset", "offset_end", "start", "end"], axis=1, inplace=True
    )
    merged_df.rename(
        columns={"new_chrom": "chrom", "new_start": "start", "new_end": "end"},
        inplace=True,
    )
    return merged_df


def load_bed_files(
    bed_dir: Path, bed_cols: List[str], split_bed: Optional[Path] = None
) -> pd.DataFrame:
    df_list = []
    for bed_i in bed_dir.glob("*.bed"):
        logger.info(f"Load {bed_i} ...")
        sample_name = bed_i.stem.rstrip(".cov")
        bed_i_columns = [*bed_cols, sample_name]
        df_i = pd.read_table(bed_i, header=None, names=[sample_name], usecols=[3])
        df_list.append(df_i)
    df = reduce(
        lambda x, y: pd.merge(x, y, left_index=True, right_index=True),
        df_list,
    )
    if split_bed is None:
        return df
    else:
        return merge_chr(df, split_bed)


def main(
    cds_cov_dir: Path,
    out_file: Path,
    transcript_id: bool = False,
    span: int = 120,
) -> None:
    if transcript_id:
        loci_columns = BED_COLUMNS[:]
    else:
        loci_columns = BED_COLUMNS[:-1]
    cds_df = load_bed_files(cds_cov_dir, loci_columns)
    df_list = []
    for cov in READS_COV:
        df_matrix_bool = cds_df >= (cov * span)
        cover_df = df_matrix_bool.sum()
        cover_ratio_df = cover_df / df_matrix_bool.shape[0]
        cover_ratio_df.name = f"coverage_{cov}x"
        df_list.append(cover_ratio_df)
    merged_df = pd.concat(df_list, axis=1)
    merged_df.to_csv(out_file)


if __name__ == "__main__":
    typer.run(main)
