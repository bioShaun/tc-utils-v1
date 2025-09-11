from pathlib import Path

import pandas as pd
import typer
from pyfaidx import Fasta
from tqdm import tqdm


def extract_fasta_id(fasta_file: Path) -> pd.DataFrame:
    fasta_id_list = []
    for record in tqdm(Fasta(fasta_file)):
        fasta_id_list.append({"id": record.name})
    return pd.DataFrame(fasta_id_list)


def main(
    fasta_file: Path,
    minimap_dir: Path,
    out_file: Path,
    max_match_count: int = 1,
    output_all: bool = False,
) -> None:
    minimap_files = sorted(minimap_dir.glob("./*"))
    fasta_id_df = extract_fasta_id(fasta_file)
    id_count_df_list = []
    for n, minimap_file in enumerate(tqdm(minimap_files)):
        minimap_df = pd.read_table(
            minimap_file,
            usecols=[0, 9],
            names=["id", "match_len"],
        )

        id_count_df = minimap_df["id"].value_counts()
        if not output_all:
            id_count_df = id_count_df[id_count_df <= max_match_count]
            minimap_df = minimap_df[minimap_df["id"].isin(id_count_df.index)]
        id_count_df = id_count_df.reset_index()
        id_count_df.columns = ["id", "minimap_match"]

        best_match_df = minimap_df.groupby("id")["match_len"].max().reset_index()
        best_match_df.columns = ["id", "best_match_len"]
        id_count_df = pd.merge(id_count_df, best_match_df, on="id", how="left")
        id_count_df["best_match_len"] = (
            id_count_df["best_match_len"].fillna(0).astype(int)
        )
        id_count_df_list.append(id_count_df)

    id_count_df = pd.concat(id_count_df_list)
    id_count_df = pd.merge(id_count_df, fasta_id_df, on="id", how="left")
    id_count_df["minimap_match"].fillna(0, inplace=True)
    id_count_df["minimap_match"] = id_count_df["minimap_match"].astype(int)
    id_count_df["best_match_len"].fillna(0, inplace=True)
    id_count_df["best_match_len"] = id_count_df["best_match_len"].astype(int)
    id_count_df.to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)
