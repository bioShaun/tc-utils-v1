from pathlib import Path
from typing import Optional

import pandas as pd
import typer


def main(bed: Path, split_bed: Path, restore_bed: Optional[Path] = None) -> None:
    if restore_bed is None:
        restore_bed = bed.with_suffix(".restore.bed")

    split_bed_df = pd.read_table(
        split_bed,
        header=None,
        names=["old_chrom", "split_start", "split_end", "chrom"],
    )
    bed_df = pd.read_table(bed, header=None, names=["chrom", "start", "end"])
    restore_bed_df = bed_df.merge(split_bed_df, on="chrom")
    restore_bed_df["old_start"] = (
        restore_bed_df["start"] + restore_bed_df["split_start"]
    )
    restore_bed_df["old_end"] = restore_bed_df["end"] + restore_bed_df["split_start"]
    restore_bed_df.to_csv(
        restore_bed,
        header=False,
        columns=["old_chrom", "old_start", "old_end"],
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    typer.run(main)
