#!/usr/bin/env python3
"""
Calculate aggregated effective alignment lengths (XM and MD) per read from a SAM/BAM file.

For each read, it calculates the minimum and maximum effective length across all its alignments,
and the total number of alignments. Effective length is defined as max(eff_xm, eff_md).
"""

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

import pysam
import typer
from loguru import logger
from rich.console import Console
from tqdm import tqdm

console = Console()
app = typer.Typer(help="聚合计算 SAM/BAM 文件每条 read 的有效比对长度（取 XM 和 MD 的最大值）")


@dataclass
class ReadStats:
    """Stores aggregated statistics for a single read."""

    min_eff: int | None = None
    max_eff: int | None = None
    count: int = 0

    def update(self, eff_len: int) -> None:
        """Update min/max and increment alignment count."""
        if self.min_eff is None:
            self.min_eff = self.max_eff = eff_len
        else:
            self.min_eff = min(self.min_eff, eff_len)
            self.max_eff = max(self.max_eff, eff_len)
        self.count += 1


def validate_file(path: Path, name: str) -> None:
    """Validate that a file exists and is a file."""
    if not path.exists():
        logger.error(f"{name} file not found: {path}")
        console.print(f"[red]Error:[/red] {name} file not found: {path}")
        raise typer.Exit(code=1)

    if not path.is_file():
        logger.error(f"Path is not a file: {path}")
        console.print(f"[red]Error:[/red] Path is not a file: {path}")
        raise typer.Exit(code=1)


def parse_md_mismatches(md_tag: str) -> int:
    """解析 MD:Z 标签，返回 mismatch 个数（排除删除 ^ 符号）"""
    # MD 标签中数字代表匹配，字母代表错配，^ 表示删除
    # 例如 10A5^AC2 代表 10个匹配，1个错配A，5个匹配，删除AC，2个匹配
    # 我们只关心错配（字母），不计入删除
    mismatches = re.findall(r"(?<!\^)[A-Z]", md_tag)
    return len(mismatches)


def compute_alignment_eff_len(read: pysam.AlignedSegment) -> int | None:
    """
    计算单条比对记录的有效长度：max(eff_xm, eff_md)
    eff_xm = CIGAR M 长度 - XM
    eff_md = CIGAR M 长度 - MD mismatches
    """
    if read.is_unmapped or read.cigartuples is None:
        return None

    # op 0 is 'M' (alignment match/mismatch)
    cigar_m_len = sum(length for (op, length) in read.cigartuples if op == 0)

    eff_xm = None
    if read.has_tag("XM"):
        eff_xm = cigar_m_len - read.get_tag("XM")

    eff_md = None
    if read.has_tag("MD"):
        md_mis = parse_md_mismatches(read.get_tag("MD"))
        eff_md = cigar_m_len - md_mis

    # Return the maximum of the two valid effective lengths
    valid_lengths = [v for v in (eff_xm, eff_md) if v is not None]
    return max(valid_lengths) if valid_lengths else None


def process_sam(sam_path: Path) -> dict[str, ReadStats]:
    """遍历 SAM 文件并按 readID 聚合统计信息"""
    stats_map: dict[str, ReadStats] = {}

    logger.info(f"Opening SAM/BAM file: {sam_path}")
    with pysam.AlignmentFile(str(sam_path), "r") as sam:
        # fetch(until_eof=True) ensures we see all records including unmapped
        for read in tqdm(sam.fetch(until_eof=True), desc="Processing alignments", unit="alns"):
            qname = read.query_name
            if qname not in stats_map:
                stats_map[qname] = ReadStats()

            eff_len = compute_alignment_eff_len(read)
            if eff_len is not None:
                stats_map[qname].update(eff_len)
            else:
                # Still count unmapped records if we want to track total reads,
                # but the update logic above only counts valid alignments.
                # If we want align_count to be 'valid alignments', we keep it as is.
                # If we want align_count to be 'total records', we move self.count += 1.
                # Based on requirement "unmapped report count=0", we keep count for mapped only.
                pass

    return stats_map


def write_stats(stats_map: dict[str, ReadStats], out_path: Path) -> None:
    """将聚合结果写入 TSV 文件"""
    logger.info(f"Writing results to: {out_path}")
    with out_path.open("w") as f:
        f.write("readID\tmin_eff_len\tmax_eff_len\talign_count\n")
        for qname, stats in tqdm(stats_map.items(), desc="Writing output", unit="reads"):
            min_str = str(stats.min_eff) if stats.min_eff is not None else "NA"
            max_str = str(stats.max_eff) if stats.max_eff is not None else "NA"
            f.write(f"{qname}\t{min_str}\t{max_str}\t{stats.count}\n")


@app.command()
def main(
    sam_file: Annotated[Path, typer.Argument(help="输入 SAM/BAM 文件路径")],
    output: Annotated[
        Path | None,
        typer.Option("--output", "-o", help="输出 TSV 文件路径 (默认: <input>.eff_len.tsv)"),
    ] = None,
    verbose: Annotated[bool, typer.Option("--verbose", "-v", help="启用详细日志")] = False,
) -> None:
    """
    聚合计算 SAM/BAM 每条 read 的有效长度统计。

    算法：
    1. 有效长度 = max(CIGAR_M - XM, CIGAR_M - MD_mismatches)
    2. 按 readID 聚合，输出最小/最大有效长度及有效比对次数。
    """
    # 配置日志
    log_level = "DEBUG" if verbose else "INFO"
    logger.remove()
    logger.add(
        lambda msg: console.print(msg, end=""),
        level=log_level,
        format="<level>{message}</level>",
    )

    # 验证输入
    validate_file(sam_file, "SAM/BAM")

    # 确定输出路径
    output_path = output if output else sam_file.with_suffix(".eff_len.tsv")

    # 处理数据
    try:
        stats_map = process_sam(sam_file)
        write_stats(stats_map, output_path)
        logger.success(f"✓ 处理完成！结果已写入: {output_path}")
    except Exception as e:
        logger.exception("处理过程中发生错误")
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
