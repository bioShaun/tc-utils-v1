import pandas as pd
import typer
from pathlib import Path


def main(gt: Path, compare: Path, out: Path, dp: Path, min_depth: int = 10) -> None:
    gt_df = pd.read_table(gt)
    dp_df = pd.read_table(dp, header=None, names=list(gt_df.columns))
    dp_df.replace(".", 0, inplace=True)
    out_list = []
    consist_df_list = []
    in_consist_df_list = []
    compare_df = pd.read_table(compare, header=None, names=["comp1", "comp2"])
    for i in compare_df.itertuples():
        compare_name = f"{i.comp1}|{i.comp2}"
        compare_df = gt_df[["CHROM", "POS", i.comp1, i.comp2]].copy()
        compare_dp_df = dp_df[[i.comp1, i.comp2]].astype("int")
        rm_low_dp_df = compare_dp_df[
            (compare_dp_df[i.comp1] >= min_depth)
            & (compare_dp_df[i.comp2] >= min_depth)
        ]
        compare_df = compare_df.loc[rm_low_dp_df.index]
        rm_na_df = compare_df[
            (compare_df[i.comp1] != "./.") & (compare_df[i.comp1] != "./.")
        ].copy()
        consistent_count = len(rm_na_df[rm_na_df[i.comp1] == rm_na_df[i.comp2]])
        consist_df_list.append(
            rm_na_df[rm_na_df[i.comp1] == rm_na_df[i.comp2]][["CHROM", "POS"]].copy()
        )
        in_consist_df_list.append(
            rm_na_df[rm_na_df[i.comp1] != rm_na_df[i.comp2]][["CHROM", "POS"]].copy()
        )
        consistent_rate = consistent_count / len(rm_na_df)
        inconsistent_count = len(rm_na_df) - consistent_count
        inconsistent_rate = inconsistent_count / len(rm_na_df)
        out_list.append(
            {
                "compare": compare_name,
                "repeatability_count": consistent_count,
                "repeatability_rate": consistent_rate,
            }
        )
    consist_df = pd.concat(consist_df_list)
    consist_df = pd.DataFrame(consist_df.groupby(["CHROM", "POS"]).size())
    consist_df.columns = ["consist_count"]
    in_consist_df = pd.concat(in_consist_df_list)
    in_consist_df = pd.DataFrame(in_consist_df.groupby(["CHROM", "POS"]).size())
    in_consist_df.columns = ["in_consist_count"]
    site_df = consist_df.merge(in_consist_df, left_index=True, right_index=True)
    site_df.fillna(0, inplace=True)
    site_df = site_df.astype("int")
    df = pd.DataFrame(out_list)
    df.to_excel(out, index=False)
    site_df.to_csv(out.with_suffix(".sites.tsv"), sep="\t")


if __name__ == "__main__":
    typer.run(main)
