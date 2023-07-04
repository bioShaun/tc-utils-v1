from functools import partial
from operator import attrgetter
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
import typer
from loguru import logger
from pandarallel import pandarallel
from tqdm import tqdm


def get_score(
    row, gc_bias: float, gc_score: int, mapping_score: int, extreme_gc_limit: float
) -> float:
    score = 10 * (4 - row["score"])
    score_low = 0.5 - gc_bias
    score_high = 0.5 + gc_bias
    if row["GC"] >= score_low and row["GC"] <= score_high:
        score += gc_score
    else:
        if row["GC"] < 0.5:
            if row["GC"] < score_low - extreme_gc_limit:
                add_score = (
                    -1
                    * (((score_low - row["GC"] - extreme_gc_limit) / 0.1) ** 2)
                    * gc_score
                )
            else:
                add_score = (1 - ((row["GC"] - score_low) / 0.1) ** 2) * gc_score
        else:
            if row["GC"] > score_high + extreme_gc_limit:
                add_score = (
                    -1
                    * (((row["GC"] - score_high - extreme_gc_limit) / 0.1) ** 2)
                    * gc_score
                )
            else:
                add_score = (1 - ((score_high - row["GC"]) ** 2)) * gc_score
        score += add_score
    allele_count = len(str(row["allele"]).split(","))
    isIndel = len(str(row["allele"])) > 3
    if isIndel:
        score -= 50
    if allele_count < 2:
        allele_count = 2
    if "maf" in row:
        if row["maf"] < 0.1:
            score -= 100
        if row["maf"] < 0.2:
            score -= 50
    mapping_score_0 = mapping_score * abs(row["match_0"] - 1)
    mapping_score_60 = mapping_score * abs(row["match_60"] - 1) * 1.5
    mapping_score_100 = mapping_score * abs(row["match_100"] - 1) ** 2
    mapping_score = -1 * (mapping_score_0 + mapping_score_60 + mapping_score_100)
    allele_count_score = -1 * (allele_count - 2) * 5
    try:
        return int(score + mapping_score + allele_count_score)
    except:
        print(row)
        print(score, mapping_score, allele_count_score)
        raise ValueError("score error")


def make_window(chr_size_df: pd.DataFrame, region_size: int) -> pd.DataFrame:
    df_list = []
    for row in chr_size_df.itertuples():
        for i in range(0, row.chrom_length, region_size):
            df_list.append([row.Index, i, i + region_size])
    return pd.DataFrame(df_list, columns=["chrom", "start", "end"])


def choose_best_row(rows):
    return sorted(rows, key=attrgetter("score2"), reverse=True)[0]


def filter_by_dist(df: pd.DataFrame, snp_min_dist: int) -> pd.DataFrame:
    passed_rows = []
    candidate_rows = []
    for row in df.itertuples():
        if len(candidate_rows) == 0:
            candidate_rows.append(row)
            continue
        if (
            row.chrom == candidate_rows[0].chrom
            and row.pos - candidate_rows[-1].pos < snp_min_dist
        ):
            candidate_rows.append(row)
            continue
        else:
            passed_rows.append(choose_best_row(candidate_rows))
            candidate_rows = [row]
    if len(candidate_rows) > 0:
        passed_rows.append(choose_best_row(candidate_rows))

    return pd.DataFrame(passed_rows).drop(columns=["Index"])


def select_by_region(
    chr_size: Path,
    region_size: int,
    output_prefix: Path,
    candidate_files: Optional[List[str]] = typer.Option(None),
    region_count: int = 1,
    high_score: int = 0,
    auto_add: bool = False,
    threads: int = 8,
    gc_bias: float = 0.15,
    gc_score: int = 1000,
    mapping_score: int = 500,
    extreme_gc_limit: float = 0.1,
    snp_min_dist: int = 50,
    total_select: int = typer.Option(None, help="total select number"),
    score_cutoff: int = 0,
) -> None:
    # reg_df = pd.read_csv(
    #     region_cover, sep="\t", header=None, names=["region", "chrom", "pos"]
    # )
    # reg_df = reg_df[reg_df.chrom != "."]
    if candidate_files is None:
        raise ValueError("candidate_files is required")
    pandarallel.initialize(nb_workers=threads)
    logger.info("loading candidates...")
    can_df = pd.concat(
        [pd.read_csv(candidate_file) for candidate_file in candidate_files]
    )
    can_df["chrom"] = can_df["chrom"].astype(str)
    chr_df = pd.read_table(
        chr_size, header=None, names=["chrom", "chrom_length"], index_col=0
    )
    chr_df.index = chr_df.index.astype(str)
    window_df = make_window(chr_df, region_size)
    add_region_df_list = []
    for chrom, each_chrom_df in can_df.groupby("chrom"):
        logger.info(f"add region for chrom {chrom}")
        if chrom not in chr_df.index:
            continue
        chrom_len = chr_df.loc[chrom, "chrom_length"]  # type: ignore
        each_chrom_cut = pd.cut(
            each_chrom_df["pos"], range(0, chrom_len + region_size, region_size)
        )
        each_chrom_df["region"] = each_chrom_cut.map(
            lambda x: f"{chrom}_{x.left}_{x.right}"
        )
        add_region_df_list.append(each_chrom_df)
    add_region_df = pd.concat(add_region_df_list)

    logger.info("calculate score...")
    my_get_score = partial(
        get_score,
        gc_bias=gc_bias,
        gc_score=gc_score,
        mapping_score=mapping_score,
        extreme_gc_limit=extreme_gc_limit,
    )
    # add_region_df["score2"] = add_region_df.apply(my_get_score, axis=1)
    add_region_df["score2"] = add_region_df.parallel_apply(my_get_score, axis=1)  # type: ignore
    add_region_df.sort_values("score2", ascending=False, inplace=True)
    add_region_df.drop_duplicates(subset=["chrom", "pos"], inplace=True)
    center_df = add_region_df[add_region_df.sequence_type == "center"].copy()
    region_dfs = []

    logger.info("select by chrom: ")
    # for row in tqdm(window_df.itertuples(), total=len(window_df)):
    for region, df in tqdm(center_df.groupby("region")):
        df = filter_by_dist(df, snp_min_dist)
        high_score_df = df[df["score"] <= high_score]
        low_score_df = df[df["score"] > high_score]
        if len(high_score_df) >= region_count:
            region_dfs.append(high_score_df)
        else:
            region_dfs.append(high_score_df)
            region_dfs.append(low_score_df[: region_count - len(high_score_df)])

    final_df = pd.concat(region_dfs)
    final_df.sort_values(["chrom", "pos"], inplace=True)
    # final_df = final_df[final_df["score2"] >= score_cutoff]
    # final_df = filter_by_dist(final_df, snp_min_dist)
    if auto_add and total_select is not None and len(final_df) < total_select:
        missed_count = total_select - len(final_df)
        filter_out_df = center_df[~center_df["id"].isin(final_df["id"])].copy()
        tmp_df = filter_out_df[:missed_count]
        minimal_score = tmp_df["score2"].min()
        add_df = filter_out_df[filter_out_df["score2"] >= minimal_score].copy()
        add_df = add_df.sample(n=missed_count, random_state=1)
        auto_add_file = f"{output_prefix}.auto_add.csv"
        add_df.to_csv(auto_add_file, index=False, float_format="%.3f")
        final_df = pd.concat([final_df, add_df])
    candidate_file = f"{output_prefix}.select.csv"
    region_count_file = f"{output_prefix}.region_count.tsv"

    region_count_df = pd.DataFrame(final_df.groupby("region").size())
    region_count_df.columns = ["count"]
    region_count_df.reset_index(inplace=True)
    region_count_df[["chrom", "start", "end"]] = region_count_df["region"].str.split(
        "_", expand=True
    )
    region_count_df[["start", "end"]] = region_count_df[["start", "end"]].astype(int)
    region_count_df = region_count_df[["chrom", "start", "end", "count"]]

    region_count_df = window_df.merge(region_count_df, how="left").fillna(0)
    region_count_df["count"] = region_count_df["count"].astype("int")
    region_count_df.to_csv(region_count_file, sep="\t", index=False, header=False)
    final_df.to_csv(candidate_file, index=False, float_format="%.3f")


if __name__ == "__main__":
    typer.run(select_by_region)
