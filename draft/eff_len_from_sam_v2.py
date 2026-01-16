#!/usr/bin/env python3
"""
Calculate aggregated effective alignment lengths (XM and MD) per read from a SAM/BAM file.
Optimized for large files using name-sorting and streaming processing to minimize memory usage.

For each read, it calculates the minimum and maximum effective length across all its alignments,
and the total number of alignments. Effective length is defined as max(eff_xm, eff_md).
"""

import os
import re
import tempfile
from dataclasses import dataclass
from itertools import groupby
from pathlib import Path
from typing import Annotated, Iterable

import pysam
import typer
from loguru import logger
from rich.console import Console
from tqdm import tqdm

console = Console()
app = typer.Typer(
    help="聚合计算 SAM/BAM 文件每条 read 的有效比对长度（取 XM 和 MD 的最大值），针对大文件进行 O(1) 内存优化。"
)


@dataclass(slots=True)
class ReadStats:
    """Stores aggregated statistics for a single read."""

    min_eff: int | None = None
    max_eff: int | None = None
    count: int = 0

    def update(self, eff_len: int) -> None:
        """Update min/max and increment alignment count."""
        if self.min_eff is None or self.max_eff is None:
            self.min_eff = self.max_eff = eff_len
        else:
            if eff_len < self.min_eff:
                self.min_eff = eff_len
            if eff_len > self.max_eff:
                self.max_eff = eff_len
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
        tag_val = read.get_tag("XM")
        if isinstance(tag_val, (int, float)):
            eff_xm = cigar_m_len - int(tag_val)

    eff_md = None
    if read.has_tag("MD"):
        tag_val = read.get_tag("MD")
        if isinstance(tag_val, str):
            md_mis = parse_md_mismatches(tag_val)
            eff_md = cigar_m_len - md_mis

    # Return the maximum of the two valid effective lengths
    valid_lengths: list[int] = [v for v in (eff_xm, eff_md) if v is not None]
    return max(valid_lengths) if valid_lengths else None


def detect_sort_order(sam_path: Path) -> str:
    """检测 SAM/BAM 排序方式"""
    try:
        with pysam.AlignmentFile(str(sam_path), "r") as sam:
            header = sam.header.to_dict()
            return header.get("HD", {}).get("SO", "unknown")
    except Exception as e:
        logger.warning(f"Could not detect sort order: {e}")
        return "unknown"


def aggregate_group(alignments: Iterable[pysam.AlignedSegment]) -> ReadStats:
    """聚合一组相同 readID 的比对记录"""
    stats = ReadStats()
    for read in alignments:
        eff_len = compute_alignment_eff_len(read)
        if eff_len is not None:
            stats.update(eff_len)
    return stats


def process_streaming(sam_path: Path, out_path: Path, threads: int = 1) -> None:
    """流式处理已按名称排序的 SAM/BAM 文件"""
    logger.info(f"Processing alignments from: {sam_path}")

    with pysam.AlignmentFile(str(sam_path), "r", threads=threads) as sam:
        with out_path.open("w") as f:
            f.write("readID\tmin_eff_len\tmax_eff_len\talign_count\n")

            # 使用 groupby 进行流式聚合，前提是文件已按 query_name 排序
            # fetch(until_eof=True) 确保读取所有记录，包括未比对的
            it = groupby(sam.fetch(until_eof=True), key=lambda r: r.query_name)

            for qname, group in tqdm(it, desc="Streaming processing", unit="reads"):
                stats = aggregate_group(group)
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
    threads: Annotated[int, typer.Option("--threads", "-t", help="排序和读取使用的线程数")] = 4,
    keep_sorted: Annotated[bool, typer.Option("--keep-sorted", help="保留生成的临时排序文件")] = False,
    verbose: Annotated[bool, typer.Option("--verbose", "-v", help="启用详细日志")] = False,
) -> None:
    """
    聚合计算 SAM/BAM 每条 read 的有效长度统计。
    自动检测排序状态，如果未按名称排序，将自动调用 pysam 内部排序。
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

    # 检测排序
    sort_order = detect_sort_order(sam_file)
    logger.info(f"Detected sort order: {sort_order}")

    temp_sorted_bam = None
    processing_file = sam_file

    if sort_order != "queryname":
        logger.info("File is not name-sorted. Sorting now (pysam sort -n)...")
        # 创建临时排序文件
        temp_dir = output_path.parent
        temp_fd, temp_path = tempfile.mkstemp(suffix=".name_sorted.bam", dir=temp_dir)
        os.close(temp_fd)
        temp_sorted_bam = Path(temp_path)

        try:
            pysam.sort("-n", "-@", str(threads), "-o", str(temp_sorted_bam), str(sam_file))
            logger.info(f"Sorting complete: {temp_sorted_bam}")
            processing_file = temp_sorted_bam
        except Exception as e:
            logger.error(f"Sorting failed: {e}")
            if temp_sorted_bam.exists():
                temp_sorted_bam.unlink()
            raise typer.Exit(code=1)

    # 流式处理
    try:
        process_streaming(processing_file, output_path, threads=threads)
        logger.success(f"✓ 处理完成！结果已写入: {output_path}")
    except Exception as e:
        logger.exception("处理过程中发生错误")
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(code=1)
    finally:
        # 清理临时文件
        if temp_sorted_bam and temp_sorted_bam.exists() and not keep_sorted:
            logger.info(f"Removing temporary sorted file: {temp_sorted_bam}")
            temp_sorted_bam.unlink()


if __name__ == "__main__":
    app()
