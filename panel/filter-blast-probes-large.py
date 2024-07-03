from pathlib import Path

import pandas as pd
import typer
from tqdm import tqdm


def main(
    blast_dir: Path,
    out_file: Path,
    min_match_length: int = 24,
    max_match_count: int = 1,
) -> None:
    blast_files = sorted(blast_dir.glob("./*"))
    result_df = pd.DataFrame()

    for n, blast_file in enumerate(tqdm(blast_files)):
        blast_df = pd.read_table(
            blast_file, usecols=[0, 2, 3], names=["id", "identity", "match_len"]
        )
        blast_df = blast_df[blast_df["match_len"] >= min_match_length]
        id_count_df = blast_df["id"].value_counts()
        filter_id_df = id_count_df[id_count_df <= max_match_count].reset_index()
        filter_id_df.columns = ["id", "blast_match"]
        if n == 0:
            result_df = filter_id_df
        else:
            result_df = result_df.append(filter_id_df, ignore_index=True)

    result_df.to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)
