import typer

import pandas as pd

from loguru import logger
from tqdm import tqdm
from pathlib import Path
from pandarallel import pandarallel

from typing import List, Optional


def get_score(row) -> float:
    score = 0
    if row["GC"] >= 0.35 and row["GC"] <= 0.65:
        score += 100
    else:
        if row["GC"] < 0.5:
            punish_score = (0.35 - row["GC"]) * 100
        else:
            punish_score = (row["GC"] - 0.65) * 100
        score -= punish_score
    allele_count = len(str(row["allele"]).split(","))
    isIndel = len(str(row["allele"])) > 3
    if isIndel:
        score -= 5
    if allele_count < 2:
        allele_count = 2
    if "maf" in row:
        if row["maf"] < 0.1:
            score -= 100
        if row["maf"] < 0.2:
            score -= 50
    return (
        score
        - 10 * (row["match_0"] - 1)
        + 1000 * (4 - row["score"])
        - (allele_count - 2) * 0.1
    )


def make_window(chr_size_df: pd.DataFrame, region_size: int) -> pd.DataFrame:
    df_list = []
    for row in chr_size_df.itertuples():
        for i in range(0, row.chrom_length, region_size):
            df_list.append([row.Index, i, i + region_size])
    return pd.DataFrame(df_list, columns=["chrom", "start", "end"])


def select_by_region(
    chr_size: Path,
    region_size: int,
    output_prefix: Path,
    candidate_files: Optional[List[str]] = typer.Option(None),
    region_count: int = 1,
    high_score: int = 0,
    auto_add: bool = False,
    threads: int = 8,
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
    chr_df = pd.read_table(
        chr_size, header=None, names=["chrom", "chrom_length"], index_col=0
    )
    window_df = make_window(chr_df, region_size)
    add_region_df_list = []
    for chrom, each_chrom_df in can_df.groupby("chrom"):
        logger.info(f"add region for chrom {chrom}")
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
    add_region_df["score2"] = add_region_df.parallel_apply(get_score, axis=1)  # type: ignore
    add_region_df.sort_values("score2", ascending=False, inplace=True)
    add_region_df.drop_duplicates(subset=["chrom", "pos"], inplace=True)
    center_df = add_region_df[add_region_df.sequence_type == "center"].copy()
    region_dfs = []

    missed_count = 0
    region_count_df_list = []

    logger.info("select by chrom: ")
    # for row in tqdm(window_df.itertuples(), total=len(window_df)):
    for region, df in tqdm(center_df.groupby("region")):
        # df = center_df[center_df["region"] == row.region]
        high_score_df = df[df["score"] <= high_score]
        low_score_df = df[df["score"] > high_score]
        if len(df) < region_count:
            missed_count += region_count - len(df)
        if len(high_score_df) >= region_count:
            region_dfs.append(high_score_df)
            real_region_count = len(high_score_df)
        else:
            region_dfs.append(high_score_df)
            region_dfs.append(low_score_df[: region_count - len(high_score_df)])
            real_region_count = len(high_score_df) + len(
                low_score_df[: region_count - len(high_score_df)]
            )
        chrom, start, end = region.split("_")  # type: ignore
        region_count_df_list.append([chrom, int(start), int(end), real_region_count])
    final_df = pd.concat(region_dfs)
    if auto_add:
        filter_out_df = center_df[~center_df["id"].isin(final_df["id"])].copy()
        tmp_df = filter_out_df[:missed_count]
        minimal_score = tmp_df["score2"].min()
        add_df = filter_out_df[filter_out_df["score2"] >= minimal_score].copy()
        add_df = add_df.sample(n=missed_count, random_state=1)
        auto_add_file = output_prefix.with_suffix(".auto_add.csv")
        add_df.to_csv(auto_add_file, index=False, float_format="%.3f")
        final_df = pd.concat([final_df, add_df])
    final_df.sort_values(["chrom", "pos"], inplace=True)
    candidate_file = f"{output_prefix}.select.csv"
    region_count_file = f"{output_prefix}.region_count.tsv"
    region_count_df = pd.DataFrame(
        region_count_df_list, columns=["chrom", "start", "end", "count"]
    )
    region_count_df = window_df.merge(region_count_df, how="left").fillna(0)
    region_count_df["count"] = region_count_df["count"].astype("int")
    region_count_df.to_csv(region_count_file, sep="\t", index=False, header=False)
    final_df.to_csv(candidate_file, index=False, float_format="%.3f")


if __name__ == "__main__":
    typer.run(select_by_region)
