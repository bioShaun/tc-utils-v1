import typer
import pandas as pd

from pathlib import Path

RAW_HEADER = ["chrom", "pos", "allele", "score"]
BED_HEADER = [
    "chrom",
    "pos",
]


def isSnp(alleles):
    return len(alleles) == 3


def main(
    candidates_csv: Path, bed: Path, out_file: Path, score: int = 1, snp: bool = True
) -> None:
    candidates_df = pd.read_csv(candidates_csv, header=None, names=RAW_HEADER)
    bed_df = pd.read_csv(bed, header=None, names=BED_HEADER, usecols=[0, 2], sep="\t")
    candidates_df = candidates_df.merge(bed_df)
    candidates_df["score"] = score
    if snp:
        candidates_df["is_snp"] = candidates_df["allele"].map(isSnp)
        candidates_df = candidates_df[candidates_df["is_snp"]]
        candidates_df.drop("is_snp", axis=1, inplace=True)
    candidates_df.to_csv(out_file, index=False, header=False)


if __name__ == "__main__":
    typer.run(main)
