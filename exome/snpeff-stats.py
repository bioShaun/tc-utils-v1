from pathlib import Path
from typing import List, Tuple
import typer
import pandas as pd
from io import StringIO


def load_stats(snpeff_file: Path, stats="effects") -> pd.DataFrame:
    flag = False
    string_list = []
    with snpeff_file.open() as snpeff_inf:
        for eachline in snpeff_inf:
            if f"Count by {stats}" in eachline:
                flag = True
                continue
            if flag:
                if eachline.strip():
                    if eachline.startswith("#"):
                        flag = False
                    else:
                        string_list.append(eachline)
    table_content = "".join(string_list)
    df = pd.read_csv(StringIO(table_content), sep=" , ")
    return df


def reformat_df(
    df: pd.DataFrame, sample_name: str
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    sample_df = df.copy()
    sample_df["sample_id"] = sample_name
    # sample_df["stats"] = sample_df.apply(
    #     lambda x: f'{x["Count"]} ({x["Percent"]})', axis=1
    # )
    count_df = (
        sample_df[["Type", "sample_id", "Count"]]
        .set_index(["Type", "sample_id"])
        .unstack(0)
    )
    count_df.columns = count_df.columns.droplevel().values
    count_df["Total"] = df["Count"].sum()
    percent_df = (
        sample_df[["Type", "sample_id", "Percent"]]
        .set_index(["Type", "sample_id"])
        .unstack(0)
    )
    percent_df.columns = percent_df.columns.droplevel().values
    return (count_df, percent_df)  # type: ignore


def main(
    all_stats: Path, sample_stats: Path, out_path: Path, stats: str = "effects"
) -> None:
    all_stats_df = load_stats(all_stats, stats)
    sample_stats_files = sample_stats.glob("*.stat.csv")
    sample_count_df_list = []
    sample_percent_df_list = []
    for each_file in sample_stats_files:
        sample_name = each_file.name.rstrip(".stat.csv")
        sample_df = load_stats(each_file, stats)
        sample_count_df, sample_percent_df = reformat_df(sample_df, sample_name)
        sample_count_df_list.append(sample_count_df)
        sample_percent_df_list.append(sample_percent_df)
    merged_sample_count_df = pd.concat(sample_count_df_list)
    merged_sample_percent_df = pd.concat(sample_percent_df_list)
    merged_sample_count_df.fillna(0, inplace=True)
    merged_sample_percent_df.fillna("0%", inplace=True)
    count_columns = ["Total", *merged_sample_percent_df.columns]

    with pd.ExcelWriter(out_path) as writer:
        all_stats_df.to_excel(writer, sheet_name="ALL", index=False)
        merged_sample_count_df.to_excel(
            writer, sheet_name="Sample Count", columns=count_columns
        )
        merged_sample_percent_df.to_excel(writer, sheet_name="Sample Percent")


if __name__ == "__main__":
    typer.run(main)
