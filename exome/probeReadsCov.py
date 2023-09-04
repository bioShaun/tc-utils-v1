from pathlib import Path
import typer
import pandas as pd
from functools import reduce


def main(cov_dir: Path, out_file: Path) -> None:
    df_list = []
    for file_i in cov_dir.glob("*.txt"):
        sample_name = file_i.name.replace(".cov.txt", "")
        df = pd.read_table(
            file_i, header=None, names=["chrom", "start", "end", "read_count"]
        )
        df["probe_length"] = df["end"] - df["start"]
        df[sample_name] = df["read_count"] / df["probe_length"]
        df_list.append(df[["chrom", "start", "end", sample_name]].copy())
    merged_df = reduce(
        lambda x, y: pd.merge(x, y, on=["chrom", "start", "end"]),
        df_list,
    )
    merged_df.to_csv(out_file, sep="\t", index=False, float_format="%.1f")


if __name__ == "__main__":
    typer.run(main)
