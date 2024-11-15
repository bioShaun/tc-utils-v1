from pathlib import Path
from typing import Optional

import pandas as pd
import typer
from tqdm import tqdm


def main(
    blast_dir: Path,
    out_file: Path,
    min_match_length: int = 24,
    max_match_count: int = 1,
    output_all: bool = False,
    show_max_length: bool = False,
    min_identity: Optional[int] = None,
    max_mismatch_count: Optional[int] = None,
    max_gap_count: Optional[int] = None,
) -> None:
    blast_files = sorted(blast_dir.glob("./*"))

    for n, blast_file in enumerate(tqdm(blast_files)):
        blast_df = pd.read_table(
            blast_file,
            usecols=[0, 2, 3, 4, 5],
            names=["id", "identity", "match_len", "mismatch", "gapopen"],
        )
        if not output_all:
            blast_df = blast_df[blast_df["match_len"] >= min_match_length]
            if min_identity is not None:
                blast_df = blast_df[blast_df["identity"] >= min_identity]
            if max_mismatch_count is not None:
                blast_df = blast_df[blast_df["mismatch"] <= max_mismatch_count]
            if max_gap_count is not None:
                blast_df = blast_df[blast_df["gapopen"] <= max_gap_count]
        id_count_df = blast_df["id"].value_counts()
        if not output_all:
            id_count_df = id_count_df[id_count_df <= max_match_count]
        id_count_df = id_count_df.reset_index()
        id_count_df.columns = ["id", "blast_match"]
        if show_max_length:
            id_max_length_df = blast_df.groupby("id")["match_len"].max()
            id_max_length_df = id_max_length_df.reset_index()
            id_max_length_df.columns = ["id", "max_length"]
            id_count_df = pd.merge(id_count_df, id_max_length_df, on="id")
        mode = "w" if n == 0 else "a"
        header = n == 0
        id_count_df.to_csv(out_file, sep="\t", index=False, mode=mode, header=header)


if __name__ == "__main__":
    typer.run(main)
