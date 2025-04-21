from functools import partial
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

    my_get_real_match_length = partial(get_real_match_length, probe_lenth=probe_length)

    for n, blast_file in enumerate(tqdm(blast_files)):
        blast_df = pd.read_table(
            blast_file,
            usecols=[0, 2, 4, 5, 6, 7],
            names=[
                "id",
                "identity",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
            ],
        )

        blast_df["match_len"] = blast_df["qend"] - blast_df["qstart"] + 1
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

        blast_df = blast_df.groupby("id").head(2)
        blast_df["real_match_length"] = (
            blast_df["match_len"] - blast_df["mismatch"] - blast_df["gapopen"]
        )
        best_match_idx = blast_df.groupby("id")["real_match_length"].idxmax()
        best_match_df = blast_df[blast_df.index.isin(best_match_idx)]
        best_match_gap_df = best_match_df.groupby("id")["gapopen"].max().reset_index()
        best_match_gap_df.columns = ["id", "best_mathc_gapopen"]
        best_match_mismatch_df = (
            best_match_df.groupby("id")["mismatch"].max().reset_index()
        )
        best_match_mismatch_df.columns = ["id", "best_match_mismatch"]
        id_count_df = pd.merge(id_count_df, best_match_gap_df, on="id", how="left")
        id_count_df = pd.merge(id_count_df, best_match_mismatch_df, on="id", how="left")
        other_match_df = blast_df[~blast_df.index.isin(best_match_idx)]
        second_max_length_df = (
            other_match_df.groupby("id")["real_match_length"].max().reset_index()
        )
        second_max_length_df.columns = ["id", "second_max_length"]
        second_best_gap_df = other_match_df.groupby("id")["gapopen"].max().reset_index()
        second_best_gap_df.columns = ["id", "second_best_gapopen"]
        second_best_mismatch_df = (
            other_match_df.groupby("id")["mismatch"].max().reset_index()
        )
        second_best_mismatch_df.columns = ["id", "second_best_mismatch"]
        id_count_df = pd.merge(id_count_df, second_max_length_df, on="id", how="left")
        id_count_df = pd.merge(id_count_df, second_best_gap_df, on="id", how="left")
        id_count_df = pd.merge(
            id_count_df, second_best_mismatch_df, on="id", how="left"
        )
        id_count_df["second_max_length"] = (
            id_count_df["second_max_length"].fillna(0).astype(int)
        )
        id_count_df["second_best_gapopen"] = (
            id_count_df["second_best_gapopen"].fillna(0).astype(int)
        )
        id_count_df["second_best_mismatch"] = (
            id_count_df["second_best_mismatch"].fillna(0).astype(int)
        )

        if show_max_length:
            id_max_length_df = (
                blast_df.groupby("id")["real_match_length"].max().reset_index()
            )
            id_max_length_df.columns = ["id", "max_length"]
            id_count_df = pd.merge(id_count_df, id_max_length_df, on="id")
        mode = "w" if n == 0 else "a"
        header = n == 0
        id_count_df.to_csv(out_file, sep="\t", index=False, mode=mode, header=header)


if __name__ == "__main__":
    typer.run(main)
