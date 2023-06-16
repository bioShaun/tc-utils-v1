from functools import partial
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
import typer
from Bio import SeqIO

CHUNK_SIZE = 1_000_000


def add_seq(row: pd.Series, record_dict: Dict[str, Any], half_length: int):
    pos = row["pos"] - 1
    start = pos - half_length
    if start < 0:
        start = 0
    end = pos + 1 + half_length
    chrom = str(row["chrom"])
    left_seq = record_dict[chrom].seq[start:pos]
    right_seq = record_dict[chrom].seq[pos + 1 : end]
    middle_seq = f"[{row['ref']}/{row['alt']}]"
    return f"{left_seq}{middle_seq}{right_seq}"


def main(vcf: Path, ref: Path, out_dir: Path, half_length: int = 200) -> None:
    ref_dict = {}
    for record in SeqIO.parse(ref, "fasta"):
        ref_dict[record.id] = record
    dfs = pd.read_table(
        vcf,
        chunksize=CHUNK_SIZE,
        header=None,
        usecols=[0, 1, 3, 4],
        names=["chrom", "pos", "ref", "alt"],
        comment="#",
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    seq_table = out_dir / "seq.table.csv"
    seq_fasta = out_dir / "seq.fasta"
    seq_fasta_inf = open(seq_fasta, "a")
    my_add_seq = partial(add_seq, record_dict=ref_dict, half_length=half_length)

    for n, df in enumerate(dfs):
        mask = (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)
        df = df[mask].copy()
        df["sequence"] = df.apply(my_add_seq, axis=1)
        header = n == 0
        df.to_csv(seq_table, header=header, index=False, mode="a")
        for row in df.itertuples():
            seq_id = f"{row.chrom}_{row.pos}"
            seq_fasta_inf.write(f">{seq_id}\n{row.sequence}\n")


if __name__ == "__main__":
    typer.run(main)
