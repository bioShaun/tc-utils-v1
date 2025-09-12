from functools import partial
from pathlib import Path

import pandas as pd
import typer
from pyfaidx import Fasta
from tqdm import tqdm

tqdm.pandas(desc="my bar!")


def get_flank_seq(
    genome_fa: Fasta, flank_size: int, chrom: str, pos: int, ref: str, alt: str
):
    left = genome_fa[chrom][pos - 1 - flank_size : pos - 1]
    right = genome_fa[chrom][pos : pos + flank_size]
    return f"{left}[{ref}/{alt}]{right}"


def main(
    design_table: Path,
    genome_path: Path,
    flank_size: int = 100,
):
    df = pd.read_table(design_table)
    out_file = design_table.with_suffix(".flank_seq.tsv")
    genome_fa = Fasta(genome_path, as_raw=True)
    my_get_flank_seq = partial(
        get_flank_seq, genome_fa=genome_fa, flank_size=flank_size
    )
    df["flank_seq"] = df.progress_apply(
        lambda x: my_get_flank_seq(
            chrom=x["chrom"], pos=x["pos"], ref=x["ref"], alt=x["alt"]
        ),
        axis=1,
    )
    df.to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)
