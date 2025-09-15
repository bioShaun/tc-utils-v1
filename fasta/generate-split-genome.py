import re
from pathlib import Path

import pandas as pd
import pyranges as pr
import scipy as sp
import typer
from numpy import split

SPLIT_SIZE = 500_000_000


def repeat_not_in_gene(
    repeat_df: pd.DataFrame, feature_df: pd.DataFrame
) -> pd.DataFrame:
    gr1 = pr.PyRanges(repeat_df)
    gr2 = pr.PyRanges(feature_df)
    overlaped_df = gr1.join(gr2)
    filter_repeart_df = repeat_df[~repeat_df["id"].isin(overlaped_df.df["id"])]
    return filter_repeart_df


def generate_split_chr_bed(best_split_site_df: pd.DataFrame) -> pd.DataFrame:
    split_df_list = []
    for row in best_split_site_df.itertuples():
        chrom = row.Chromosome
        start = row.Start
        end = row.End
        split_point = (start + end) // 2
        split_df_list.append(
            {
                "Chromosome": chrom,
                "start": 0,
                "end": split_point,
                "id": f"{chrom}a",
            }
        )
        split_df_list.append(
            {
                "Chromosome": chrom,
                "start": split_point,
                "end": row.chrom_size,
                "id": f"{chrom}b",
            }
        )
    return pd.DataFrame(split_df_list)


def split_chrom(
    chr_df: pd.DataFrame,
    repeat_df: pd.DataFrame,
    repeat_min_size: int = 100,
) -> pd.DataFrame:
    add_repeat_df = repeat_df.merge(chr_df)
    add_repeat_df["split_start"] = add_repeat_df["chrom_size"] * 0.3
    add_repeat_df["split_end"] = add_repeat_df["chrom_size"] * 0.7
    split_site_add_repeat_df = add_repeat_df[
        (add_repeat_df["Start"] >= add_repeat_df["split_start"])
        & (add_repeat_df["End"] <= add_repeat_df["split_end"])
        & (add_repeat_df["repeat_size"] >= repeat_min_size)
    ]
    best_split_site_idx = split_site_add_repeat_df.groupby("Chromosome")[
        "repeat_size"
    ].idxmax()
    best_split_site_df = split_site_add_repeat_df.loc[best_split_site_idx].copy()
    # check if there is any chromosome not in best_split_site_df
    not_in_best_split_site = chr_df[
        ~chr_df["Chromosome"].isin(best_split_site_df["Chromosome"])
    ]
    if not_in_best_split_site.empty:
        print("All chromosomes have split sites.")
        return generate_split_chr_bed(best_split_site_df)
    raise ValueError(
        f"There are {len(not_in_best_split_site)} chromosomes: {not_in_best_split_site['Chromosome'].to_list()} not in best_split_site_df."
    )


def main(genome_fai: Path, repeat_out: Path, gff: Path, split_cat_bed: Path) -> None:
    fai_df = pd.read_table(
        genome_fai,
        header=None,
        names=["Chromosome", "chrom_size"],
        usecols=[0, 1],
        sep="\t",
    )
    need_to_split_df = fai_df[fai_df["chrom_size"] >= SPLIT_SIZE]
    do_not_need_to_split_df = fai_df[fai_df["chrom_size"] < SPLIT_SIZE].copy()

    repeat_df = pd.read_table(
        repeat_out,
        header=None,
        sep=r"\s+",
        usecols=[4, 5, 6, 12],
        names=["Chromosome", "Start", "End", "repeat_size"],
    )
    repeat_df["id"] = (
        repeat_df["Chromosome"].astype(str)
        + "_"
        + repeat_df["Start"].astype(str)
        + "_"
        + repeat_df["End"].astype(str)
    )

    feature_df = pd.read_table(
        gff,
        header=None,
        sep="\t",
        names=["Chromosome", "Start", "End"],
        usecols=[0, 3, 4],
    )

    filtered_repeat_df = repeat_not_in_gene(repeat_df, feature_df)
    split_chr_bed = split_chrom(need_to_split_df, filtered_repeat_df)
    do_not_need_to_split_df["start"] = 0
    do_not_need_to_split_df["end"] = do_not_need_to_split_df["chrom_size"]
    do_not_need_to_split_df["id"] = do_not_need_to_split_df["Chromosome"]
    split_chr_bed_all = pd.concat([split_chr_bed, do_not_need_to_split_df])
    split_chr_bed_all.to_csv(
        split_cat_bed,
        sep="\t",
        index=False,
        header=False,
        columns=split_chr_bed.columns,
    )


if __name__ == "__main__":
    typer.run(main)
