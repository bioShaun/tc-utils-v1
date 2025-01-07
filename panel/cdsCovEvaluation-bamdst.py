from functools import reduce
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import typer
from loguru import logger

BED_COLUMNS = ["chrom", "start", "end", "transcript_id"]


def merge_chr(df: pd.DataFrame, split_bed: Path) -> pd.DataFrame:
    split_bed_df = pd.read_csv(
        split_bed,
        header=None,
        names=["new_chrom", "offset", "offset_end", "chrom"],
        sep="\t",
    )
    print(df)
    print(split_bed_df)
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
    bed_dir: Path, cov_cutoff: Optional[float] = None
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    df_list = []
    bed_list = sorted(list(bed_dir.glob("*/region.tsv.gz")))
    bed_df = pd.read_table(bed_list[0], usecols=[0, 1, 2])
    for bed_i in bed_list:
        logger.info(f"Load {bed_i} ...")
        sample_name = bed_i.parent.name
        df_i = pd.read_table(bed_i, usecols=[3])
        df_i.columns = [sample_name]
        if cov_cutoff is not None:
            if df_i[sample_name].quantile() < cov_cutoff:
                continue
        df_list.append(df_i)
    df = reduce(
        lambda x, y: pd.merge(x, y, left_index=True, right_index=True),
        df_list,
    )
    return bed_df, df


def get_stats_df(df: pd.DataFrame) -> pd.DataFrame:
    min_cov = df.min(1)
    max_cov = df.max(1)
    ave_cov = df.sum(1) / df.shape[1]
    quantile_cov = df.quantile(0.5, axis=1)
    merged_df = pd.concat([min_cov, max_cov, ave_cov, quantile_cov], axis=1)
    merged_df.columns = ["min_cov", "max_cov", "mean_cov", "quantile_cov"]
    return merged_df


def main(
    cds_cov_dir: Path,
    out_file: Path,
    cov: List[int] = [1, 5, 10, 20, 30, 50, 100],
    split_bed: Path = typer.Option(None),
    cov_cutoff: float = typer.Option(None),
) -> None:
    bed_df, df_matrix = load_bed_files(cds_cov_dir, cov_cutoff=cov_cutoff)
    stats_df = get_stats_df(df_matrix)
    cov_df_list = []
    for cov_i in cov:
        cov_i_df_matrix = df_matrix >= cov_i
        print("coverd samples")
        cover_df = cov_i_df_matrix.sum(1)
        print("cover ratio")
        cover_ratio_df = cover_df / cov_i_df_matrix.shape[1]
        cover_ratio_df.name = f"coverage_{cov_i}x"
        cov_df_list.append(cover_ratio_df)
    if not split_bed is None:
        bed_df = merge_chr(bed_df, split_bed)
    cover_ratio_df = pd.concat([bed_df, stats_df, *cov_df_list], axis=1)
    cover_ratio_df.to_csv(out_file, index=False, float_format="%.3f", sep="\t")


if __name__ == "__main__":
    typer.run(main)
