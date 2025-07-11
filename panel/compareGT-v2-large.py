from enum import Enum
from itertools import combinations
from pathlib import Path
from typing import Optional

import delegator
import pandas as pd
import typer
from loguru import logger
from tqdm import tqdm

LOCATION_COLS = ["CHROM", "POS", "REF", "ALT"]


class InputType(str, Enum):
    VCF = "vcf"
    TABLE = "table"


def vcf2gt(vcf_file: Path, force: bool = False) -> Path:
    gt_file = vcf_file.with_suffix(".gt.txt.gz")
    if gt_file.exists() and not force:
        return gt_file
    gt_file.parent.mkdir(parents=True, exist_ok=True)
    cmd = f'bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n" {vcf_file} | sed -re "s;\\|;/;g" | gzip > {gt_file}'
    logger.info(f"run: {cmd}")
    delegator.run(cmd)
    return gt_file


def get_sample_names(vcf_file: Path) -> list:
    cmd = f"bcftools query -l {vcf_file}"
    logger.info(f"run: {cmd}")
    return delegator.run(cmd).out.strip().split("\n")


def map_gt(gt: str) -> str:
    """Map a genotype string to its corresponding state."""
    if gt == "./.":
        return "MISS"
    ref, alt = gt.split("/")
    if ref == alt:
        return "HOM"
    return "HET"


def compare_gt(df: pd.DataFrame, compare_a: str, compare_b: str):
    stats_gt_df = df[[compare_a, compare_b]].map(map_gt)
    mask1 = stats_gt_df[compare_a] != "MISS"
    mask2 = stats_gt_df[compare_b] != "MISS"
    non_miss_df = stats_gt_df[mask1 & mask2]
    non_miss_count = len(non_miss_df)
    a_hom_count = len(non_miss_df[non_miss_df[compare_a] == "HOM"])
    a_het_count = len(non_miss_df[non_miss_df[compare_a] == "HET"])
    b_hom_count = len(non_miss_df[non_miss_df[compare_b] == "HOM"])
    b_het_count = len(non_miss_df[non_miss_df[compare_b] == "HET"])

    mask3 = non_miss_df[compare_a] == "HOM"
    mask4 = non_miss_df[compare_b] == "HOM"
    homo_non_miss_df = non_miss_df[mask3 & mask4]
    het_non_miss_df = non_miss_df[~non_miss_df.index.isin(homo_non_miss_df.index)]
    total_homo_count = len(homo_non_miss_df)
    total_het_count = len(het_non_miss_df)
    homo_compare_df = df[df.index.isin(homo_non_miss_df.index)]
    het_compare_df = df[df.index.isin(het_non_miss_df.index)]
    homo_equal_count = len(
        homo_compare_df[homo_compare_df[compare_a] == homo_compare_df[compare_b]]
    )
    homo_equal_percent = (
        0 if total_homo_count == 0 else 100 * homo_equal_count / total_homo_count
    )

    het_equal_count = len(
        het_compare_df[het_compare_df[compare_a] == het_compare_df[compare_b]]
    )
    het_equal_percent = (
        0 if total_het_count == 0 else 100 * het_equal_count / total_het_count
    )
    total_a_equal_b = homo_equal_count + het_equal_count
    total_a_equal_b_percent = (
        0 if non_miss_count == 0 else 100 * total_a_equal_b / non_miss_count
    )
    total_a_not_equal_b = non_miss_count - total_a_equal_b
    total_a_not_equal_b_percent = (
        0 if non_miss_count == 0 else 100 * total_a_not_equal_b / non_miss_count
    )
    return [
        compare_a,
        compare_b,
        len(df),
        non_miss_count,
        a_hom_count,
        a_het_count,
        b_hom_count,
        b_het_count,
        total_a_equal_b,
        round(total_a_equal_b_percent, 3),
        homo_equal_count,
        round(homo_equal_percent, 3),
        het_equal_count,
        round(het_equal_percent, 3),
        total_a_not_equal_b,
        round(total_a_not_equal_b_percent, 3),
    ]


def main(
    input_file: Path,
    output_file: Path,
    force: bool = False,
    input_type: InputType = InputType.VCF,
    sample_file: Optional[Path] = typer.Option(None),
    compare_list: Optional[Path] = typer.Option(None),
):
    if input_type == InputType.VCF:
        gt_file = vcf2gt(input_file, force=force)
        sample_list = get_sample_names(input_file)
    else:
        gt_file = input_file
        if sample_file is None:
            raise ValueError("Must provide sample_list if input_type is not VCF")
        sample_list = pd.read_csv(sample_file, header=None)[0].tolist()
    df = pd.read_table(gt_file, header=None, names=[*LOCATION_COLS, *sample_list])
    if compare_list is not None:
        compare_df = pd.read_table(compare_list, header=None, names=["A", "B"])
    else:
        compare_dict = [{"A": a, "B": b} for a, b in combinations(sample_list, 2)]
        compare_df = pd.DataFrame(compare_dict)

    with open(output_file, "w", encoding="utf-8") as out_inf:
        out_inf.write(
            "A,B,总位点数,有效位点,A_纯合,A_杂合,B_纯合,B_杂合,整体相似度,整体相似度%,纯合相似度,纯合相似度%,杂合相似度,杂合相似度%,整体差异位点数,整体差异%\n"
        )
        for row in tqdm(compare_df.itertuples(), total=len(compare_df)):
            row_a, row_b = str(row.A), str(row.B)
            if row_a in df.columns and row_b in df.columns:
                # sample_stats_list.append(compare_gt(df, row.A, row.B, human=human))  # type: ignore
                out_line_list = [str(each) for each in compare_gt(df, row_a, row_b)]
                out_line_str = ",".join(out_line_list)
                out_inf.write(f"{out_line_str}\n")


if __name__ == "__main__":
    typer.run(main)
