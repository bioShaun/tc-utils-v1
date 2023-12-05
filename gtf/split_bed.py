from pathlib import Path
import pandas as pd
import typer


def main(bed: Path, split_bed: Path) -> None:
    bed_df = pd.read_table(bed, header=None, names=["chrom", "start", "end"])
    split_bed_df = pd.read_table(
        split_bed, header=None, names=["chrom", "split_start", "split_end", "new_chrom"]
    )
    merge_df = bed_df.merge(split_bed_df)
    merge_df = merge_df[
        (merge_df["start"] >= merge_df["split_start"])
        & (merge_df["start"] < merge_df["split_end"])
    ]
    merge_df["new_start"] = merge_df["start"] - merge_df["split_start"]
    merge_df["new_end"] = merge_df["end"] - merge_df["split_start"]
    out_bed = bed.with_suffix(".split.bed")
    merge_df.to_csv(
        out_bed, header=False, columns=["new_chrom", "new_start", "new_end"]
    )


if __name__ == "__main__":
    typer.run(main)
