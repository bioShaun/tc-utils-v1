#!/usr/bin/env python3
import re
from typing import Optional, Tuple

import pysam
import typer
from tqdm import tqdm

app = typer.Typer(help="计算 Bowtie2 SAM 文件每条 read 的有效比对长度（XM 和 MD）")

# ---------------- Core Functions ---------------- #


def parse_md_mismatches(md_tag: str) -> int:
    """解析MD:Z标签，返回 mismatch 个数"""
    mismatches = re.findall(r"[A-Z]", md_tag)
    return len(mismatches)


def compute_effective_lengths(read) -> Tuple[Optional[int], Optional[int]]:
    """
    返回一条 read 的有效比对长度：
    - eff_xm: CIGAR M长度 - XM（None 如果没有 XM）
    - eff_md: CIGAR M长度 - MD mismatch（None 如果没有 MD）
    """
    # unmapped 或 CIGAR 缺失的情况
    if read.cigartuples is None:
        return None, None

    cigar_len = sum(length for (op, length) in read.cigartuples if op == 0)

    eff_xm = None
    if read.has_tag("XM"):
        eff_xm = cigar_len - read.get_tag("XM")

    eff_md = None
    if read.has_tag("MD"):
        md_mis = parse_md_mismatches(read.get_tag("MD"))
        eff_md = cigar_len - md_mis

    return eff_xm, eff_md


def process_sam(sam_path: str, out_path: str):
    """逐条读取 SAM 文件并计算 XM 和 MD 有效长度，直接写入文件"""
    with pysam.AlignmentFile(sam_path, "r") as sam, open(out_path, "w") as f:
        f.write("readID\teff_len_XM\teff_len_MD\n")
        for read in tqdm(sam.fetch(until_eof=True), desc="Processing SAM"):
            eff_xm, eff_md = compute_effective_lengths(read)
            xm_str = str(eff_xm) if eff_xm is not None else "NA"
            md_str = str(eff_md) if eff_md is not None else "NA"
            f.write(f"{read.query_name}\t{xm_str}\t{md_str}\n")


# ---------------- CLI ---------------- #


@app.command()
def main(sam_file: str, out_file: str = "effective_length.tsv"):
    """
    sam_file: 输入SAM文件
    out_file: 输出TSV文件，格式：readID \t eff_len_XM \t eff_len_MD
    """
    process_sam(sam_file, out_file)
    typer.echo(f"Done! 结果已写入 {out_file}")


if __name__ == "__main__":
    app()
    app()
