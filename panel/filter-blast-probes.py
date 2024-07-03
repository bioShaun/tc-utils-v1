from pathlib import Path

import pandas as pd
import typer


def main(
    blast_dir: Path, out_file: Path, min_len: int = 25, max_match_count: int = 1
) -> None:
    df_list = [
        pd.read_table(each, usecols=[0, 2, 3], names=["id", "identity", "match_len"])
        for each in blast_dir.glob("./*")
    ]
    merged_df = pd.concat(df_list)
    merged_df = merged_df[merged_df["match_len"] >= min_len]
    id_count_df = merged_df["id"].value_counts()
    filter_id_df = id_count_df[id_count_df <= max_match_count].reset_index()
    filter_id_df.columns = ["id", "blast_match"]
    filter_id_df.to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)
