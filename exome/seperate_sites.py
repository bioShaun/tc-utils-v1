from pathlib import Path

import pandas as pd
import typer
from loguru import logger


def get_indel_type(row: pd.Series) -> str:
    ref, alt = row["ref"], row["alt"]
    if len(ref) > len(alt):
        return "DEL"
    if len(ref) < len(alt):
        return "INS"
    return "SNP"


def get_index_len(row: pd.Series) -> int:
    ref, alt = row["ref"], row["alt"]
    ref_len = len(ref)
    alt_list = [len(x) for x in alt.split(",")]
    alt_len_max, alt_len_min = max(alt_list), min(alt_list)
    return max(abs(ref_len - alt_len_max), abs(ref_len - alt_len_min))


def indel_right_pos(row: pd.Series) -> int:
    if row["variant_type"] in ["DEL", "MIXED"]:
        return row["pos"] + row["indel_length"]
    return row["pos"]


def main(mixed_table: Path, output_dir: Path):
    df = pd.read_csv(mixed_table, sep="\t")
    df["ref"] = df["alleles"].str.split(",").str[0]
    df["alt"] = df["alleles"].str.split(",").str[1]
    df["variant_type"] = df.apply(get_indel_type, axis=1)
    df["indel_length"] = df.apply(get_index_len, axis=1)
    df["id"] = df.apply(lambda x: f"{x['chrom']}_{x['pos']}", axis=1)
    snp_df = df[df["variant_type"] == "SNP"]
    indel_df = df[df["variant_type"] != "SNP"]
    snp_df.drop(columns=["ref", "alt", "indel_length"]).to_csv(
        output_dir / "snp.tsv", sep="\t", index=False
    )
    indel_df.drop(columns=["ref", "alt", "indel_length"]).to_csv(
        output_dir / "indel.left.tsv", sep="\t", index=False
    )
    indel_right_df = indel_df.copy()
    indel_right_df["pos"] = indel_right_df.apply(indel_right_pos, axis=1)
    indel_right_df.drop(columns=["ref", "alt", "indel_length"]).to_csv(
        output_dir / "indel.right.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    typer.run(main)
