from typing import List, Optional
import typer

import pandas as pd
from enum import Enum
from pathlib import Path


class OutFmt(str, Enum):
    GT = "GT"
    AB = "AB"


def parentGT(row: pd.Series, cols: List[str]) -> Optional[str]:
    rowValues = list(set([row[col] for col in cols]))
    if len(rowValues) > 1 or rowValues[0] in ["./.", "0/1"]:
        return None
    else:
        return rowValues[0]


def convertGT(row: pd.Series, **kwargs) -> str:
    child_name = kwargs.pop("child_name", None)
    miss_fmt = kwargs.pop("miss_fmt", None)
    out_fmt = kwargs.pop("out_fmt", None)
    if row[child_name] == "./.":
        return miss_fmt
    if row[child_name] == row["P1"]:
        return "AA"
    elif row[child_name] == row["P2"]:
        return "BB"
    else:
        return "AB"


def gt2linkagemap(
    gt_file: Path,
    out_file: Path,
    p1: str,
    p2: str,
    out_fmt=typer.Option("AB", help="AB/GT"),
    miss_fmt=typer.Option("N", help="N/NN/--"),
) -> None:
    gt_df = pd.read_csv(gt_file, sep="\t")
    p1_names = p1.split(",")
    p2_names = p2.split(",")
    gt_df["P1"] = gt_df.apply(parentGT, axis=1, args=(p1_names,))
    gt_df["P2"] = gt_df.apply(parentGT, axis=1, args=(p2_names,))
    gt_df.dropna(inplace=True)
    gt_df = gt_df[gt_df["P1"] != gt_df["P2"]]
    child_names = [
        each
        for each in gt_df.columns[4:]
        if each not in [*p1_names, *p2_names, "P1", "P2"]
    ]
    for child in child_names:
        gt_df[child] = gt_df.apply(
            convertGT, axis=1, child_name=child, miss_fmt=miss_fmt, out_fmt=out_fmt
        )
    gt_df["P1"] = "AA"
    gt_df["P2"] = "BB"
    for p1 in p1_names:
        gt_df[p1] = "AA"
    for p2 in p2_names:
        gt_df[p2] = "BB"

    loci_columns = gt_df.columns[:4]
    out_cols = [*loci_columns, *p1_names, *p2_names, *child_names]
    gt_df.to_csv(out_file, sep="\t", index=False, columns=out_cols)

    child_stats_count = (
        gt_df[child_names].applymap(str).apply(pd.Series.value_counts).T.fillna(0)
    )
    child_stats_ratio = (child_stats_count.T / child_stats_count.sum(1)).T
    child_stats_count_df = child_stats_count.reset_index()
    child_stats_count_df["stats"] = "Count"
    child_stats_ratio_df = child_stats_ratio.reset_index()
    child_stats_ratio_df["stats"] = "Ratio"
    sample_stats_all = pd.concat([child_stats_count_df, child_stats_ratio_df])
    merged_sample_stats = (
        sample_stats_all.set_index(["index", "stats"])
        .melt(ignore_index=False)
        .reset_index()
        .sort_values(["stats", "variable"])
        .set_index(["index", "variable", "stats"])
        .unstack([2, 1])
    )
    merged_sample_stats.columns = merged_sample_stats.columns.droplevel()
    merged_sample_stats["Count"] = merged_sample_stats["Count"].fillna(0).astype("int")
    merged_sample_stats.index.name = ""
    merged_sample_stats.columns.names = ["", ""]
    child_stats_file = out_file.with_name(f"gt.bysample.{out_file.name}")
    merged_sample_stats.to_csv(child_stats_file, sep="\t", float_format="%.3f")

    by_loci_df = gt_df[[*loci_columns, *child_names]].set_index([*loci_columns])
    by_loci_count = by_loci_df.T.apply(pd.Series.value_counts).T.fillna(0)
    by_loci_ratio = (by_loci_count.T / by_loci_count.sum(1)).T
    by_loci_count_df = by_loci_count.reset_index()
    by_loci_count_df["stats"] = "Count"
    by_loci_ratio_df = by_loci_ratio.reset_index()
    by_loci_ratio_df["stats"] = "Ratio"
    loci_stats_all = pd.concat([by_loci_count_df, by_loci_ratio_df])
    merged_loci_stats = (
        loci_stats_all.set_index([*loci_columns, "stats"])
        .melt(ignore_index=False)
        .reset_index()
        .sort_values(["stats", "variable"])
        .set_index([*loci_columns, "variable", "stats"])
        .unstack([5, 4])
    )
    merged_loci_stats.columns = merged_loci_stats.columns.droplevel()
    merged_loci_stats["Count"] = merged_loci_stats["Count"].astype("int")
    merged_loci_stats.index.name = ""
    merged_loci_stats.columns.names = ["", ""]
    merged_loci_stats_file = out_file.with_name(f"gt.byLoci.{out_file.name}")
    merged_loci_stats.to_csv(merged_loci_stats_file, sep="\t", float_format="%.3f")


if __name__ == "__main__":
    typer.run(gt2linkagemap)
