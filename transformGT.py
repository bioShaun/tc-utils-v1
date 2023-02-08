from pathlib import Path
import typer
import pandas as pd
from functools import partial


SAMPLE_COL_INDEX = 5
INFO_COLUMS = ["ID", "CHROM", "POS", "REF", "ALT"]

GT_MAP = {
    "0/0": "AA",
    "1/1": "BB",
    "0/1": "AB",
    "./.": "NC",
}

app = typer.Typer()


def load_df(filename: Path, col_name: str) -> pd.DataFrame:
    df = pd.read_csv(filename, sep="\t")
    return df.melt(
        id_vars=df.columns[:SAMPLE_COL_INDEX],
        value_vars=df.columns[SAMPLE_COL_INDEX:],
        value_name=col_name,
        var_name="name",
    )


def load_parent_df(
    filename: Path, col_name: str, sample_a: str, sample_b: str
) -> pd.DataFrame:
    df = pd.read_csv(filename, sep="\t")
    id_vars = [*INFO_COLUMS, sample_a, sample_b]
    value_vars = [each for each in df.columns if each not in id_vars]
    return df.melt(
        id_vars=id_vars,
        value_vars=value_vars,
        value_name=col_name,
        var_name="name",
    )


def transformGT(row: pd.Series, min_depth: int) -> str:
    if row["depth"] < min_depth:
        return "NC"
    else:
        return GT_MAP.get(row["genotype"], row["genotype"])


def transformByParent(
    row: pd.Series, min_depth: int, sample_a: str, sample_b: str
) -> str:
    if row["depth"] < min_depth:
        return "NC"
    elif row["genotype"] == "./.":
        return "NC"
    elif row["genotype"] == row[sample_a]:
        return "AA"
    elif row["genotype"] == row[sample_b]:
        return "BB"
    else:
        return "AB"


@app.command()
def transform_by_gt(gt: Path, dp: Path, out: Path, depth: int = 5) -> None:
    if depth <= 0:
        raise ValueError("depth must be greater than zero")
    gt_df = load_df(gt, "genotype")
    dp_df = load_df(dp, "depth")
    dp_df["depth"] = dp_df["depth"].replace(".", 0).map(lambda x: int(x))
    gt_dp_df = gt_df.merge(dp_df)
    transformGTByDP = partial(transformGT, min_depth=depth)
    gt_dp_df["genotype_ab"] = gt_dp_df.apply(transformGTByDP, axis=1)
    out_df = gt_dp_df.drop(["genotype", "depth"], axis=1).pivot(
        index=["ID", "CHROM", "POS", "REF", "ALT"], columns="name", values="genotype_ab"
    )
    out_df.to_csv(out, sep="\t")


@app.command()
def transform_by_parent(
    gt: Path, dp: Path, out: Path, sample_a: str, sample_b: str, depth: int = 5
) -> None:
    if depth <= 0:
        raise ValueError("depth must be greater than zero")
    gt_df = load_parent_df(gt, "genotype", sample_a, sample_b)
    dp_df = load_parent_df(dp, "depth", sample_a, sample_b)
    dp_df = dp_df[(dp_df[sample_a] >= depth) & (dp_df[sample_b] >= depth)]
    dp_df["depth"] = dp_df["depth"].replace(".", 0).map(lambda x: int(x))
    dp_df.drop([sample_a, sample_b], axis=1, inplace=True)
    gt_dp_df = gt_df.merge(dp_df)
    transformGTByDpPa = partial(
        transformByParent, min_depth=depth, sample_a=sample_a, sample_b=sample_b
    )
    gt_dp_df["genotype_ab"] = gt_dp_df.apply(transformGTByDpPa, axis=1)
    gt_dp_df[sample_a] = "AA"
    gt_dp_df[sample_b] = "BB"
    out_df = gt_dp_df.drop(["genotype", "depth"], axis=1).pivot(
        index=["ID", "CHROM", "POS", "REF", "ALT", sample_a, sample_b],
        columns="name",
        values="genotype_ab",
    )
    out_df.to_csv(out, sep="\t")


if __name__ == "__main__":
    app()
