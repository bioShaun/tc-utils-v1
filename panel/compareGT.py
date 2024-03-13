import pandas as pd
import typer
from pathlib import Path


def main(gt: Path, compare: Path, out: Path, dp: Path, min_depth: int = 10) -> None:
    gt_df = pd.read_table(gt)
    dp_df = pd.read_table(dp, header=None, names=list(gt_df.columns))
    out_list = []
    compare_df = pd.read_table(compare, header=None, names=["comp1", "comp2"])
    for i in compare_df.itertuples():
        compare_name = f"{i.comp1}|{i.comp2}"
        compare_df = gt_df[[i.comp1, i.comp2]].copy()
        rm_low_dp_df = dp_df[
            (dp_df[i.comp1] >= min_depth) & (dp_df[i.comp2] >= min_depth)
        ]
        compare_df = compare_df.loc[rm_low_dp_df.index]
        rm_na_df = compare_df[
            (compare_df[i.comp1] != "./.") & (compare_df[i.comp1] != "./.")
        ].copy()
        consistent_count = len(rm_na_df[rm_na_df[i.comp1] == rm_na_df[i.comp2]])
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
    df = pd.DataFrame(out_list)
    df.to_excel(out, index=False)


if __name__ == "__main__":
    typer.run(main)
