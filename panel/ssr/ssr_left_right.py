from pathlib import Path

import pandas as pd
import typer
from tqdm import tqdm

BASE_BLAST_DIR = Path("./top_seq_blast_top")
OUT_COLS = ["id", "chrom", "pos", "alleles", "label"]
OUTDIR = Path(".")


def get_left_right_df(
    df: pd.DataFrame, direction: str, genome: str, out_dir: Path, offset: int
):
    direct_df = df.copy()
    if direction == "left":
        direct_df["pos"] = direct_df.apply(lambda x: min(x[8], x[9]) + offset, axis=1)
    else:
        direct_df["pos"] = direct_df.apply(lambda x: max(x[8], x[9]) - offset, axis=1)
    direct_df["chrom"] = direct_df[1]
    direct_df["alleles"] = "-/-"
    direct_df["label"] = f"wheat_ssr_fq-{genome}"
    out_file = out_dir / f"ssr.{direction}.tsv"
    direct_df.to_csv(out_file, sep="\t", index=False, columns=OUT_COLS)


def main(
    genome: str = typer.Option(..., help="Genome"),
    blast_dir: Path = typer.Argument(default=BASE_BLAST_DIR, help="Blast dir"),
    out_dir: Path = typer.Option(default=OUTDIR, help="Output dir"),
    off_set_bp: int = typer.Option(default=10, help="Offset bp"),
) -> None:

    def get_best_and_most(test_df):
        a = (test_df.groupby([1, 8]).size() >= 2).reset_index()
        most_df = a[a[0]].drop(0, axis=1).merge(test_df).drop_duplicates(subset=[1, 8])
        best_df = test_df[:1]
        return pd.concat([most_df, best_df]).drop_duplicates(subset=[1, 8])

    df_list = []
    for region_file in tqdm(blast_dir.glob("*.tsv")):
        test_df = pd.read_table(region_file, header=None)
        test_df["raw_id"] = region_file.stem.replace(".blast", "")
        df_list.append(get_best_and_most(test_df))
    df = pd.concat(df_list)
    df["dup_count"] = df.groupby("raw_id").cumcount()
    df["id"] = df.apply(
        lambda row: (
            f"{row['raw_id']}-{row['dup_count'] + 1}"
            if row["dup_count"] > 0
            else row["raw_id"]
        ),
        axis=1,
    )

    for direction in ["left", "right"]:
        get_left_right_df(df, direction, genome, out_dir, off_set_bp)


if __name__ == "__main__":
    typer.run(main)
