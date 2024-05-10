from enum import Enum
from functools import partial
from pathlib import Path
from typing import List, Optional

import pandas as pd
import typer

app = typer.Typer()


class TableColumn(Enum):
    CHROM = "CHROM"
    POS = "POS"
    REF = "REF"
    ALT = "ALT"
    SAMPLE_NAME = "SAMPLE_NAME"
    GENOTYPE = "GENOTYPE"
    P1 = "P1"
    P2 = "P2"
    STATS = "stats"
    COUNT = "Count"
    RATIO = "Ratio"


class GT_VALUE(Enum):
    NA = "./."
    HET = "0/1"
    REF = "0/0"
    ALT = "1/1"


class AB_GT(Enum):
    P1 = "AA"
    P2 = "BB"
    P12 = "AB"


class OutFmt(str, Enum):
    GT = "GT"
    AB = "AB"


class MissFmt(str, Enum):
    N = "N"
    NN = "NN"
    DASH = "--"


def parentGT(row: pd.Series, cols: List[str]) -> Optional[str]:
    rowValues = list(set([row[col] for col in cols]))
    if len(rowValues) > 1 or rowValues[0] in [GT_VALUE.NA.value, GT_VALUE.HET.value]:
        return None
    else:
        return rowValues[0]


def convertGT(row: pd.Series, miss_fmt: MissFmt, out_fmt: OutFmt) -> str:
    GT_MAP = {
        GT_VALUE.ALT.value: f"{row[TableColumn.ALT.value]}{row[TableColumn.ALT.value]}",
        GT_VALUE.REF.value: f"{row[TableColumn.REF.value]}{row[TableColumn.REF.value]}",
    }
    if row[TableColumn.GENOTYPE.value] == GT_VALUE.NA.value:
        return miss_fmt.value
    if row[TableColumn.GENOTYPE.value] == row[TableColumn.P1.value]:
        if out_fmt == OutFmt.AB:
            return AB_GT.P1.value
        else:
            return GT_MAP[row[TableColumn.P1.value]]
    elif row[TableColumn.GENOTYPE.value] == row[TableColumn.P2.value]:
        if out_fmt == OutFmt.AB:
            return AB_GT.P2.value
        else:
            return GT_MAP[row[TableColumn.P2.value]]
    else:
        if out_fmt == OutFmt.AB:
            return AB_GT.P12.value
        else:
            return f"{row[TableColumn.REF.value]}{row[TableColumn.ALT.value]}"


def npConvertGT(row: pd.Series, miss_fmt: MissFmt) -> str:
    GT_MAP = {
        GT_VALUE.ALT.value: f"{row[TableColumn.ALT.value]}{row[TableColumn.ALT.value]}",
        GT_VALUE.REF.value: f"{row[TableColumn.REF.value]}{row[TableColumn.REF.value]}",
        GT_VALUE.HET.value: f"{row[TableColumn.REF.value]}{row[TableColumn.ALT.value]}",
        GT_VALUE.NA.value: miss_fmt.value,
    }
    return GT_MAP[row[TableColumn.GENOTYPE.value]]


@app.command()
def gt2linkagemap(
    gt_file: Path,
    out_file: Path,
    p1: str,
    p2: str,
    out_fmt: OutFmt = OutFmt.AB,
    miss_fmt: MissFmt = MissFmt.N,
) -> None:
    gt_df = pd.read_csv(gt_file, sep="\t")
    p1_names = p1.split(",")
    p2_names = p2.split(",")
    #  confirm parent gt and filter
    gt_df["P1"] = gt_df.apply(parentGT, axis=1, args=(p1_names,))
    gt_df["P2"] = gt_df.apply(parentGT, axis=1, args=(p2_names,))
    gt_df.dropna(inplace=True)
    gt_df = gt_df[gt_df["P1"] != gt_df["P2"]]

    # transform
    loci_columns = gt_df.columns[:4]
    melt_gt_df = gt_df.melt(
        id_vars=[*loci_columns, TableColumn.P1.value, TableColumn.P2.value],
        var_name=TableColumn.SAMPLE_NAME.value,
        value_name=TableColumn.GENOTYPE.value,
    )
    myConvertGT = partial(convertGT, miss_fmt=miss_fmt, out_fmt=out_fmt)
    melt_gt_df[TableColumn.GENOTYPE.value] = melt_gt_df.apply(myConvertGT, axis=1)

    convert_df = (
        melt_gt_df.drop([TableColumn.P1.value, TableColumn.P2.value], axis=1)
        .set_index([*loci_columns, TableColumn.SAMPLE_NAME.value])
        .unstack(4)
    )
    convert_df.columns = convert_df.columns.droplevel()

    out_cols = [
        each
        for each in gt_df.columns
        if each not in [*loci_columns, TableColumn.P1.value, TableColumn.P2.value]
    ]
    convert_df.to_csv(out_file, sep="\t", columns=out_cols)

    if out_fmt == OutFmt.AB:
        child_columns = [
            each for each in convert_df.columns if each not in [*p1_names, *p2_names]
        ]
        stats_df = convert_df[child_columns]

        child_stats_count = (
            stats_df.applymap(str).apply(pd.Series.value_counts).T.fillna(0)
        )
        child_stats_ratio = (child_stats_count.T / child_stats_count.sum(1)).T
        child_stats_count_df = child_stats_count.reset_index()
        child_stats_count_df[TableColumn.STATS.value] = TableColumn.COUNT.value
        child_stats_ratio_df = child_stats_ratio.reset_index()
        child_stats_ratio_df[TableColumn.STATS.value] = TableColumn.RATIO.value
        sample_stats_all = pd.concat([child_stats_count_df, child_stats_ratio_df])
        merged_sample_stats = (
            sample_stats_all.set_index(
                [TableColumn.SAMPLE_NAME.value, TableColumn.STATS.value]
            )
            .melt(ignore_index=False)
            .reset_index()
            .sort_values([TableColumn.STATS.value, "variable"])
            .set_index(
                [TableColumn.SAMPLE_NAME.value, "variable", TableColumn.STATS.value]
            )
            .unstack([2, 1])
        )
        merged_sample_stats.columns = merged_sample_stats.columns.droplevel()
        merged_sample_stats[TableColumn.COUNT.value] = (
            merged_sample_stats[TableColumn.COUNT.value].fillna(0).astype("int")
        )
        merged_sample_stats.index.name = ""
        merged_sample_stats.columns.names = ["", ""]
        child_stats_file = out_file.with_name(f"gt.bysample.{out_file.name}")
        merged_sample_stats.to_csv(child_stats_file, sep="\t", float_format="%.3f")

        by_loci_count = stats_df.T.apply(pd.Series.value_counts).T.fillna(0)
        by_loci_ratio = (by_loci_count.T / by_loci_count.sum(1)).T
        by_loci_count_df = by_loci_count.reset_index()
        by_loci_count_df[TableColumn.STATS.value] = TableColumn.COUNT.value
        by_loci_ratio_df = by_loci_ratio.reset_index()
        by_loci_ratio_df[TableColumn.STATS.value] = TableColumn.RATIO.value
        loci_stats_all = pd.concat([by_loci_count_df, by_loci_ratio_df])
        merged_loci_stats = (
            loci_stats_all.set_index([*loci_columns, TableColumn.STATS.value])
            .melt(ignore_index=False)
            .reset_index()
            .sort_values([TableColumn.STATS.value, "variable"])
            .set_index([*loci_columns, "variable", TableColumn.STATS.value])
            .unstack([5, 4])
        )
        merged_loci_stats.columns = merged_loci_stats.columns.droplevel()
        merged_loci_stats[TableColumn.COUNT.value] = merged_loci_stats[
            TableColumn.COUNT.value
        ].astype("int")
        merged_loci_stats.index.name = ""
        merged_loci_stats.columns.names = ["", ""]
        merged_loci_stats_file = out_file.with_name(f"gt.byLoci.{out_file.name}")
        merged_loci_stats.to_csv(merged_loci_stats_file, sep="\t", float_format="%.3f")


@app.command()
def npGt2linkagemap(
    gt_file: Path,
    out_file: Path,
    miss_fmt: MissFmt = MissFmt.N,
) -> None:
    gt_df = pd.read_csv(gt_file, sep="\t")

    # transform
    loci_columns = gt_df.columns[:4]
    melt_gt_df = gt_df.melt(
        id_vars=[*loci_columns],
        var_name=TableColumn.SAMPLE_NAME.value,
        value_name=TableColumn.GENOTYPE.value,
    )
    myConvertGT = partial(npConvertGT, miss_fmt=miss_fmt)
    melt_gt_df[TableColumn.GENOTYPE.value] = melt_gt_df.apply(myConvertGT, axis=1)

    convert_df = melt_gt_df.set_index(
        [*loci_columns, TableColumn.SAMPLE_NAME.value]
    ).unstack(4)
    convert_df.columns = convert_df.columns.droplevel()

    out_cols = [each for each in gt_df.columns if each not in [*loci_columns]]
    convert_df.to_csv(out_file, sep="\t", columns=out_cols)


if __name__ == "__main__":
    app()
