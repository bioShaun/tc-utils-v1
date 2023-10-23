from pathlib import Path
import typer
import pandas as pd


COLUMN_IDX_MAP = {"varfilter": [1, 13, 15]}


def get_gene_sequence(
    chrom: str,
    start: int,
    end: int,
    score_type: str,
    score_file: Path,
    out_dir: Path,
    genome_base: Path,
):
    use_cols = COLUMN_IDX_MAP[score_type]
    df = pd.read_csv(score_file, usecols=use_cols)
    df = df[(df["POS"] >= start) & (df["POS"] <= end)]
    gene_tr_df = df[["Transcript", "Gene"]].drop_duplicates()


if __name__ == "__main__":
    typer.run(get_gene_sequence)
