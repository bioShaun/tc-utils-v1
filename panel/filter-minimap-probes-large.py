from pathlib import Path

import pandas as pd
import typer
from tqdm import tqdm


def main(
    minimap_dir: Path,
    out_file: Path,
    max_match_count: int = 1,
    output_all: bool = False,
) -> None:
    minimap_files = sorted(minimap_dir.glob("./*"))

    for n, minimap_file in enumerate(tqdm(minimap_files)):
        minimap_df = pd.read_table(
            minimap_file,
            usecols=[0],
            names=["id"],
        )

        id_count_df = minimap_df["id"].value_counts()
        if not output_all:
            id_count_df = id_count_df[id_count_df <= max_match_count]
        id_count_df = id_count_df.reset_index()
        id_count_df.columns = ["id", "minimap_match"]
        mode = "w" if n == 0 else "a"
        header = n == 0
        id_count_df.to_csv(out_file, sep="\t", index=False, mode=mode, header=header)


if __name__ == "__main__":
    typer.run(main)
