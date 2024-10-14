from enum import Enum
from pathlib import Path

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


def compare_gt(df: pd.DataFrame, compare_a: str, compare_b: str, human: bool = True):
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
    # if human:
    #     return {
    #         "A": compare_a,
    #         "B": compare_b,
    #         "总位点数": len(df),
    #         "有效位点": non_miss_count,
    #         "A_纯合": a_hom_count,
    #         "A_杂合": a_het_count,
    #         "B_纯合": b_hom_count,
    #         "B_杂合": b_het_count,
    #         "整体相似度": total_a_equal_b_out,
    #         "纯合相似度": homo_equal_count_out,
    #         "杂合相似度": het_equal_count_out,
    #     }
    # return {
    #     "A": compare_a,
    #     "B": compare_b,
    #     "总位点数": len(df),
    #     "有效位点": non_miss_count,
    #     "A_纯合": a_hom_count,
    #     "A_杂合": a_het_count,
    #     "B_纯合": b_hom_count,
    #     "B_杂合": b_het_count,
    #     "整体相似度": total_a_equal_b,
    #     "整体相似度%": total_a_equal_b_percent,
    #     "纯合相似度": homo_equal_count,
    #     "纯合相似度%": homo_equal_percent,
    #     "杂合相似度": het_equal_count,
    #     "杂合相似度%": het_equal_percent,
    # }
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
        total_a_equal_b_percent,
        homo_equal_count,
        homo_equal_percent,
        het_equal_count,
        het_equal_percent,
    ]


def main(
    input_file: Path,
    compare_list: Path,
    output_file: Path,
    chrom_stats: bool = False,
    force: bool = False,
    human: bool = True,
    input_type: InputType = InputType.VCF,
    sample_file: Path = typer.Option(None),
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
    compare_df = pd.read_table(compare_list, header=None, names=["A", "B"])

    with open(output_file, "w") as out_inf:
        out_inf.write(
            "A\tB\t总位点数\t有效位点\tA_纯合\tA_杂合\tB_纯合\tB_杂合\t整体相似度\tt整体相似度%\t纯合相似度\t纯合相似度%\t杂合相似度\t杂合相似度%\n"
        )
        for row in tqdm(compare_df.itertuples(), total=len(compare_df)):
            if row.A in df.columns and row.B in df.columns:
                # sample_stats_list.append(compare_gt(df, row.A, row.B, human=human))  # type: ignore
                out_line_list = [
                    str(each) for each in compare_gt(df, row.A, row.B, human=human)
                ]
                out_line_str = "\t".join(out_line_list)
                out_inf.write(f"{out_line_str}\n")


if __name__ == "__main__":
    typer.run(main)
