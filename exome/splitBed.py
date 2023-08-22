from pathlib import Path
import pandas as pd
import numpy as np
import typer


OUT_COLUMNS = ["chrom", "start", "end"]


def save_current_bedrows(rows, out_dir: Path) -> None:
    start_loci = rows[0]
    end_loci = rows[-1]
    start_pos = f"{start_loci.chrom}_{start_loci.start}"
    if start_loci.chrom == end_loci.chrom:
        end_pos = f"{end_loci.end}"
    else:
        end_pos = f"{end_loci.chrom}_{end_loci.end}"
    out_file = out_dir / f"{start_pos}_{end_pos}.bed"
    df = pd.DataFrame(rows)
    df.to_csv(out_file, sep="\t", index=False, header=False, columns=OUT_COLUMNS)


def split_bed(bed_file: Path, out_dir: Path, split_number: int) -> None:
    split_out_dir = out_dir / bed_file.stem
    if split_out_dir.exists():
        raise ValueError(f"{split_out_dir} 已存在，请检查")
    split_out_dir.mkdir(parents=True)
    bed_df = pd.read_table(bed_file, header=None, names=["chrom", "start", "end"])
    bed_df["region_length"] = bed_df["end"] - bed_df["start"]
    bed_length = bed_df["region_length"].sum()
    bed_length_per_file = bed_length // split_number
    current_bed_list = []
    current_bed_size = 0
    for row in bed_df.itertuples():
        if current_bed_size > bed_length_per_file:
            save_current_bedrows(current_bed_list, split_out_dir)
            current_bed_list = []
            current_bed_size = 0
        current_bed_size += row.region_length
        current_bed_list.append(row)
    if current_bed_list:
        save_current_bedrows(current_bed_list, split_out_dir)


def get_genome_split_length(genome_length: int, split_number: int) -> int:
    raw_length = genome_length // split_number
    multiby = int(10 ** np.floor(np.log10(raw_length)))
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
        for i in range(0, row.chrom_length, genome_split_length):
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
        split_bed(bed_file=bed_fai, out_dir=out_path, split_number=split_number)
    else:
        split_fai(fai_file=bed_fai, out_dir=out_path, split_number=split_number)


if __name__ == "__main__":
    typer.run(main)
