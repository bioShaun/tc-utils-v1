from os import name
from pathlib import Path
from typing import Tuple

import pandas as pd
import typer


def merge_intervals(group: pd.DataFrame, start_col: str, end_col: str) -> pd.DataFrame:
    # 先按 start 排序
    group = group.sort_values(start_col).reset_index(drop=True)
    merged = []
    current_start, current_end = group.loc[0, start_col], group.loc[0, end_col]

    # Ensure the values are numeric for comparison
    current_end = pd.to_numeric(current_end)

    for i in range(1, len(group)):
        s, e = group.loc[i, start_col], group.loc[i, end_col]
        s = pd.to_numeric(s)
        e = pd.to_numeric(e)
        if s <= current_end:  # 有重叠
            current_end = max(current_end, e)
        else:  # 无重叠
            merged.append((current_start, current_end))
            current_start, current_end = s, e
    merged.append((current_start, current_end))

    return pd.DataFrame(merged, columns=[start_col, end_col])


def get_flank_start_end(
    probe_start: int, probe_end: int, flank_size: int, chr_size: int
) -> Tuple[int, int]:
    probe_len = probe_end - probe_start
    extend_len = (flank_size - probe_len) // 2
    flank_start = max(0, probe_start - extend_len)
    flank_end = min(chr_size, probe_end + extend_len)
    return flank_start, flank_end


def main(
    design_table: Path,
    genome_fai: Path,
    probe_id: str,
    out_path: Path,
    flank_size: int = 200,
) -> None:
    chrom_df = pd.read_table(
        genome_fai, header=None, names=["chrom", "chrom_size"], usecols=[0, 1]
    )
    chrom_df["chrom"] = chrom_df["chrom"].astype(str)
    if design_table.suffix == ".tsv":
        df = pd.read_csv(design_table, sep="\t")
    elif design_table.suffix == ".xlsx":
        df = pd.read_excel(design_table)
    else:
        raise ValueError("Please provide a valid design table")

    df["pos_id"] = df["chrom"].astype(str) + "_" + df["pos"].astype(str)
    df["chrom"] = df["chrom"].astype(str)
    df["chrom"] = pd.Categorical(
        df["chrom"], categories=chrom_df["chrom"].tolist(), ordered=True
    )
    id_file = out_path / (probe_id + ".id")
    rm_dup_df = df.drop_duplicates(subset=["pos_id"])
    rm_dup_df.sort_values(by=["chrom", "pos"], inplace=True)
    rm_dup_df.to_csv(id_file, sep="\t", index=False, header=False, columns=["pos_id"])
    rm_dup_df["pos_0"] = rm_dup_df["pos"] - 1

    target_bed_file = out_path / (probe_id + ".bed")
    rm_dup_df.to_csv(
        target_bed_file,
        sep="\t",
        index=False,
        header=False,
        columns=["chrom", "pos_0", "pos"],
    )

    add_chr_size_df = df.merge(chrom_df)
    add_chr_size_df[["flank_start", "flank_end"]] = add_chr_size_df.apply(
        lambda row: get_flank_start_end(
            row["probe_start"], row["probe_end"], flank_size, row["chrom_size"]
        ),
        axis=1,
        result_type="expand",
    )

    result = (
        add_chr_size_df.groupby("chrom", group_keys=True)
        .apply(
            lambda x: merge_intervals(x, "flank_start", "flank_end"),
            include_groups=False,
        )
        .reset_index(level=0)
        .reset_index(drop=True)
    )
    result.to_csv(
        out_path / (probe_id + ".snpcalling.bed"), sep="\t", index=False, header=False
    )


if __name__ == "__main__":
    typer.run(main)
