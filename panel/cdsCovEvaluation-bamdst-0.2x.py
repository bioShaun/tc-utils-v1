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
    bed_dir: Path,
    sample_list: Optional[List[str]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    df_list = []
    bed_list = sorted(list(bed_dir.glob("*/depth.tsv.gz")))
    bed_df = pd.read_table(bed_list[0], usecols=[0, 1])
    bed_df.columns = ["chrom", "end"]
    bed_df["start"] = bed_df["end"] - 1
    for bed_i in bed_list:
        logger.info(f"Load {bed_i} ...")
        sample_name = bed_i.parent.name
        if sample_list is not None:
            if sample_name not in sample_list:
                continue
        df_i = pd.read_table(bed_i, usecols=[2])
        df_i.columns = [sample_name]
        depth_02x: float = df_i[sample_name].mean() * 0.2
        df_i[sample_name] = df_i[sample_name] >= depth_02x
        df_list.append(df_i)
    df = reduce(
        lambda x, y: pd.merge(x, y, left_index=True, right_index=True),
        df_list,
    )
    return bed_df, df


def main(
    cds_cov_dir: Path,
    out_file: Path,
    split_bed: Path = typer.Option(None),
    sample_path: Optional[Path] = typer.Option(None),
) -> None:
    sample_list = None
    if sample_path is not None:
        sample_list: Optional[List[str]] = pd.read_csv(sample_path, header=None)[
            0
        ].to_list()
    bed_df, df_matrix = load_bed_files(cds_cov_dir, sample_list=sample_list)
    cover_df = df_matrix.sum(1)
    cover_ratio_df = cover_df / df_matrix.shape[1]
    cover_ratio_df.name = f"coverage_0.2x"
    if not split_bed is None:
        bed_df = merge_chr(bed_df, split_bed)
    cover_ratio_df = pd.concat([bed_df, cover_ratio_df], axis=1)

    cover_ratio_df.to_csv(out_file, index=False, float_format="%.3f", sep="\t")


if __name__ == "__main__":
    typer.run(main)
