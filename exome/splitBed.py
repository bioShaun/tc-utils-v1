from pathlib import Path
import pandas as pd
import numpy as np
import typer


def split_bed(bed_file: Path, out_dir: Path, split_number: int) -> None:
    split_out_dir = out_dir / bed_file.stem
    bed_df = pd.read_table(bed_file, header=None, names=["chrom", "start", "end"])
    record_number = len(bed_df)
    bed_number_per_file = record_number // split_number


def get_genome_split_length(genome_length: int, split_number: int) -> int:
    raw_length = genome_length // split_number
    multiby = np.floor(np.log10(raw_length))
    multier = raw_length // multiby
    return multier * multiby


def split_fai(fai_file: Path, out_dir: Path, split_number: int) -> None:
    split_out_dir = out_dir / "genome"
    if split_out_dir.exists():
        raise ValueError(f"{split_out_dir} 已存在，请检查")
    split_out_dir.mkdir(parents=True)
    fai_df = pd.read_table(
        fai_file, header=None, names=["chrom", "chrom_length"], usecols=[0, 1]
    )
    genome_length = fai_df["chrom_length"].sum()
    genome_split_length = get_genome_split_length(genome_length, split_number)

    for row in fai_df.itertuples():
        for i in range(row.chrom_length):
            end = i + genome_split_length
            if end > row.chrom_length:
                end = row.chrom_length
            file_path = split_out_dir / f"{row.chrom}_{i}_{end}.bed"
            with file_path.open("w") as bed_inf:
                bed_inf.write(f"{row.chrom}\t{i}\t{end}\n")


def main(
    bed_fai: Path, out_path: Path, split_number: int = 500, is_bed: bool = True
) -> None:
    if is_bed:
        pass
    else:
        split_fai(fai_file=bed_fai, out_dir=out_path, split_number=split_number)


if __name__ == "__main__":
    typer.run(main)
