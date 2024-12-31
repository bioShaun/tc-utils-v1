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
from tqdm import tqdm
from typing_extensions import Annotated

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


def load_bed_files(bed_files: List[Path]) -> pd.DataFrame:
    df_list = []
    for bed_i in bed_files:
        logger.info(f"Load {bed_i} ...")
        sample_name = bed_i.stem.rstrip(".cov")
        df_i = pd.read_table(
            bed_i, header=None, names=["start", "end", "depth"], usecols=[1, 2, 3]
        )
        df_i["span"] = df_i["end"] - df_i["start"]
        df_i[sample_name] = df_i["depth"] / df_i["span"]
        df_list.append(df_i[[sample_name]])
    df = reduce(
        lambda x, y: pd.merge(x, y, left_index=True, right_index=True),
        df_list,
    )
    return df


def main(
    panel_cov_dir: Path,
    out_file_prefix: Path,
    mapping_summary: Annotated[
        Optional[Path], typer.Option(help="fastp data summary")
    ] = None,
    sample_map: Annotated[
        Optional[Path], typer.Option(help="library-sampleid map file")
    ] = None,
    chunk_size: int = 100,
) -> None:
    bed_files = sorted(list(panel_cov_dir.glob("*.bed")))
    stats_df_list = []
    for i in tqdm(range(0, len(bed_files), chunk_size)):
        panel_cov_df = load_bed_files(bed_files[i : i + chunk_size])
        df_list = []

        for cov in READS_COV:
            df_matrix_bool = panel_cov_df >= cov
            cover_df = df_matrix_bool.sum()
            cover_ratio_df = cover_df / df_matrix_bool.shape[0]
            cover_ratio_df.name = f"coverage_{cov}x"
            df_list.append(cover_ratio_df)
        merged_df = pd.concat(df_list, axis=1)
        capture_df = pd.DataFrame(panel_cov_df.sum(), columns=["target_bases"])
        merged_df = merged_df.merge(capture_df, left_index=True, right_index=True)
        if mapping_summary is not None:
            mapping_df = pd.read_table(mapping_summary)
            mapping_df = mapping_df[
                [
                    "Name",
                    "total length",
                    "bases mapped (cigar)",
                    "insert size average",
                    "percentage of properly paired reads (%)",
                ]
            ].copy()
            mapping_df.columns = [
                "name",
                "total_bases",
                "mapped_bases",
                "insert_size",
                "properly_paired_bases_percentage",
            ]
            mapping_df["efficiency"] = (
                mapping_df["target_bases"] / mapping_df["mapped_bases"]
            )
            merged_df = mapping_df.merge(merged_df, left_on="name", right_index=True)
        else:
            merged_df.index.name = "name"
            merged_df = merged_df.reset_index()
        if sample_map is not None:
            sample_map_df = pd.read_table(
                sample_map, header=None, usecols=[0, 1], names=["LibId", "name"]
            )
            sample_map_df["name"] = sample_map_df["name"].map(
                lambda x: re.sub("[^a-zA-Z0-9_-]", "_", str(x).strip())
            )
            sample_map_df.drop_duplicates(inplace=True)
            merged_df = sample_map_df.merge(merged_df)
        stats_df_list.append(merged_df)
    all_stats_df = pd.concat(stats_df_list)
    all_stats_df.to_excel(f"{out_file_prefix}.xlsx", index=False)
    all_stats_df.to_csv(f"{out_file_prefix}.tsv", index=False, sep="\t")


if __name__ == "__main__":
    typer.run(main)
