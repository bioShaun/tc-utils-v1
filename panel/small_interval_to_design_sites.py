import sys
from pathlib import Path
from typing import Optional, TextIO, Tuple

import typer
from rich.console import Console
from rich.progress import track

# 常量定义
FAKE_ALLELES = "-/-"
HEADER = "chrom\tpos\talleles\tregion\n"

console = Console()


def parse_interval(line: str) -> Optional[Tuple[str, int, int]]:
    """
    解析BED文件的一行数据。

    Args:
        line: BED文件的一行数据

    Returns:
        包含染色体、起始位置和终止位置的元组，如果行为空则返回None

    Raises:
        ValueError: 如果行格式不正确
    """
    line = line.strip()
    if not line:
        return None

    try:
        parts = line.split("\t")
        if len(parts) < 3:
            raise ValueError(f"Invalid line format: {line}")

        chrom, start, end = parts[:3]
        return chrom, int(start), int(end)
    except ValueError as e:
        raise ValueError(f"Error parsing line: {line}") from e


def write_sites(out: TextIO, chrom: str, pos: int, region: str) -> None:
    """
    写入位点数据到输出文件。

    Args:
        out: 输出文件对象
        chrom: 染色体名
        pos: 位点位置
        region: 区域标识
    """
    out.write(f"{chrom}\t{pos}\t{FAKE_ALLELES}\t{region}\n")


def process_interval(interval: Tuple[str, int, int]) -> Tuple[str, list[int], str]:
    """
    处理单个区间，生成需要的位点信息。

    Args:
        interval: 包含染色体、起始位置和终止位置的元组

    Returns:
        包含染色体、位点列表和区域标识的元组
    """
    chrom, start, end = interval
    center = (start + end) // 2
    region = f"{chrom}:{start}-{end}"
    return chrom, [start, center, end], region


def main(
    small_interval_bed: Path = typer.Argument(
        ...,
        help="Input BED file containing small intervals",
        exists=True,
        file_okay=True,
        dir_okay=False,
    ),
    fake_sites: Path = typer.Argument(
        ...,
        help="Output file path for fake sites",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        "-f",
        help="Force overwrite output file if it exists",
    ),
) -> None:
    """
    从小区间BED文件生成设计位点。

    对每个输入区间，在起始、中心和终止位置生成设计位点。
    """
    try:
        # 验证输出文件
        if fake_sites.exists() and not force:
            console.print(
                f"[red]Error: Output file {fake_sites} already exists. Use --force to overwrite.[/red]"
            )
            sys.exit(1)

        # 确保输出目录存在
        fake_sites.parent.mkdir(parents=True, exist_ok=True)

        # 计算总行数用于进度显示
        total_lines = sum(
            1 for line in open(small_interval_bed, encoding="utf-8") if line.strip()
        )

        with open(small_interval_bed, "r", encoding="utf-8") as f, open(
            fake_sites, "w", encoding="utf-8"
        ) as out:

            # 写入头部
            out.write(HEADER)

            # 处理每一行，显示进度
            for line in track(
                f, total=total_lines, description="Processing intervals..."
            ):
                try:
                    interval = parse_interval(line)
                    if not interval:
                        continue

                    chrom, positions, region = process_interval(interval)

                    # 写入该区间的所有位点
                    for pos in positions:
                        write_sites(out, chrom, pos, region)

                except ValueError as e:
                    console.print(
                        f"[yellow]Warning: Skipping invalid line: {str(e)}[/yellow]"
                    )
                    continue

        console.print(
            f"[green]Successfully generated fake sites at {fake_sites}[/green]"
        )

    except Exception as e:
        console.print(f"[red]Error: {str(e)}[/red]")
        sys.exit(1)


if __name__ == "__main__":
    typer.run(main)
