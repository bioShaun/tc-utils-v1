import os
from dataclasses import dataclass
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import typer

OUT_COLUMNS = ["chrom", "start", "end"]


@dataclass
class BedClass:
    chrom: str
    start: int
    end: int


def save_current_bedrows(rows, out_dir: Path, pad_num: int, prefix_idx: str) -> None:
    start_loci = rows[0]
    end_loci = rows[-1]
    start_site = str(start_loci.start)
    end_site = str(end_loci.end)

    start_pos = f"{start_loci.chrom}_{start_site}"
    end_pos = str(end_loci.end)
    if start_loci.chrom != end_loci.chrom:
        end_pos = f"{end_loci.chrom}_{end_site}"
    out_file = out_dir / f"{prefix_idx}_{start_pos}_{end_pos}.bed"
    df = pd.DataFrame(rows)
    df.to_csv(out_file, sep="\t", index=False, header=False, columns=OUT_COLUMNS)


def split_bed(bed_file: Path, out_dir: Path, split_number: int) -> None:
    split_out_dir = out_dir / bed_file.stem
    if split_out_dir.exists():
        raise ValueError(f"{split_out_dir} 已存在，请检查")
    split_out_dir.mkdir(parents=True)
    if split_number == 1:
        os.symlink(bed_file.absolute(), split_out_dir / bed_file.name)
        return
    bed_df = pd.read_table(bed_file, header=None, names=["chrom", "start", "end"])
    bed_df["region_length"] = bed_df["end"] - bed_df["start"]
    bed_length = bed_df["region_length"].sum()
    bed_length_per_file = bed_length // split_number
    current_bed_list = []
    current_bed_size = 0
    max_size = bed_df["end"].max()
    pad_num = int(np.ceil(np.log10(max_size)))
    prefix_pad_num = int(np.ceil(np.log10(split_number))) + 1
    current_idx = 0
    for row in bed_df.itertuples():
        if current_bed_size > bed_length_per_file:
            current_idx += 1
            current_idx_prefix = str(current_idx).zfill(prefix_pad_num)
            save_current_bedrows(
                current_bed_list, split_out_dir, pad_num, current_idx_prefix
            )
            current_bed_list = []
            current_bed_size = 0
        current_bed_size += row.region_length
        current_bed_list.append(row)
    if current_bed_list:
        current_idx += 1
        current_idx_prefix = str(current_idx).zfill(prefix_pad_num)
        save_current_bedrows(
            current_bed_list, split_out_dir, pad_num, current_idx_prefix
        )


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
    # fai_df["chrom"] = fai_df["chrom"].astype("str")
    # fai_df.sort_values(by=["chrom"], inplace=True)
    genome_length = fai_df["chrom_length"].sum()
    genome_split_length = get_genome_split_length(genome_length, split_number)

    max_size = fai_df["chrom_length"].max()
    pad_num = int(np.ceil(np.log10(max_size)))

    prefix_pad_num = int(np.ceil(np.log10(split_number))) + 1

    row_list = []
    current_length = 0
    step = genome_split_length // 10
    current_idx = 0
    for row in fai_df.itertuples():
        for i in range(0, row.chrom_length, step):
            if current_length >= genome_split_length:
                current_idx += 1
                current_idx_prefix = str(current_idx).zfill(prefix_pad_num)
                save_current_bedrows(
                    row_list, split_out_dir, pad_num, current_idx_prefix
                )
                row_list = []
                current_length = 0
            end = i + step
            if end > row.chrom_length:
                end = row.chrom_length
            region_length = end - i
            current_length += region_length
            last_item = row_list[-1] if row_list else None
            if last_item and last_item.chrom == row.chrom:
                last_item.end = end
            else:
                row_list.append(BedClass(chrom=str(row.chrom), start=i, end=end))
    if row_list:
        current_idx += 1
        current_idx_prefix = str(current_idx).zfill(prefix_pad_num)
        save_current_bedrows(row_list, split_out_dir, pad_num, current_idx_prefix)


def main(
    bed_fai: Path, out_path: Path, split_number: int = 400, is_bed: bool = True
) -> None:
    if is_bed:
        split_bed(bed_file=bed_fai, out_dir=out_path, split_number=split_number)
    else:
        split_fai(fai_file=bed_fai, out_dir=out_path, split_number=split_number)


if __name__ == "__main__":
    typer.run(main)
