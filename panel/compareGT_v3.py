"""
比较 VCF 文件中样本间的基因型一致性。

使用 cyvcf2 读取 VCF 文件，numpy 矩阵运算进行高性能数据处理。
采用分块处理 + 向量化计算，支持超大 VCF 文件和超多样本。

使用示例:
    # 比较所有样本对
    python compareGT_v3.py input.vcf output.csv

    # 指定比较列表（tab分隔，两列：样本A 样本B）
    python compareGT_v3.py input.vcf output.csv -c compare_pairs.txt

    # 调整分块大小（默认 10000，内存不足时可减小）
    python compareGT_v3.py input.vcf output.csv --chunk-size 5000

    # 详细日志模式
    python compareGT_v3.py input.vcf output.csv -v

输出格式:
    CSV 文件，包含以下列：
    - A, B: 比较的两个样本名
    - 总位点数: VCF 中的变异位点总数
    - 有效位点: 两个样本都有基因型的位点数
    - A_纯合, A_杂合: 样本 A 的纯合/杂合位点数
    - B_纯合, B_杂合: 样本 B 的纯合/杂合位点数
    - 整体相似度, 整体相似度%: 基因型完全一致的位点数和百分比
    - 纯合相似度, 纯合相似度%: 双方都为纯合且一致的位点统计
    - 杂合相似度, 杂合相似度%: 含杂合位点中一致的统计
    - 整体差异位点数, 整体差异%: 基因型不一致的位点统计
"""

from collections.abc import Generator
from itertools import combinations
from pathlib import Path
from typing import Annotated

import numpy as np
import typer
from cyvcf2 import VCF
from loguru import logger
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TaskProgressColumn, TextColumn

console = Console()
app = typer.Typer(help="比较 VCF 文件中样本间的基因型一致性")

# gt_types 编码常量 (cyvcf2 定义)
GT_HOM_REF = 0  # 0/0
GT_HET = 1  # 0/1, 1/0
GT_UNKNOWN = 2  # ./.
GT_HOM_ALT = 3  # 1/1


def validate_file(path: Path, name: str) -> None:
    """验证文件是否存在且为普通文件。"""
    if not path.exists():
        logger.error(f"{name} 不存在: {path}")
        console.print(f"[red]错误:[/red] {name} 不存在: {path}")
        raise typer.Exit(code=1)
    if not path.is_file():
        logger.error(f"{name} 不是文件: {path}")
        console.print(f"[red]错误:[/red] {name} 不是文件: {path}")
        raise typer.Exit(code=1)


def get_vcf_samples(vcf_path: Path) -> list[str]:
    """获取 VCF 文件中的样本列表。"""
    vcf = VCF(str(vcf_path))
    samples = list(vcf.samples)
    vcf.close()
    return samples


def load_vcf_as_matrix(
    vcf_path: Path,
    chunk_size: int = 10000,
) -> Generator[np.ndarray, None, None]:
    """
    分块读取 VCF 文件，每块返回 numpy int8 矩阵。

    Args:
        vcf_path: VCF 文件路径
        chunk_size: 每块的变异数量

    Yields:
        (variants, samples) 形状的 int8 矩阵
        编码: 0=HOM_REF, 1=HET, 2=HOM_ALT, 3=UNKNOWN
    """
    vcf = VCF(str(vcf_path))
    chunk: list[np.ndarray] = []

    for variant in vcf:
        gt = variant.gt_types.copy()
        # 将 multi-allelic 位点标记为 UNKNOWN，排除出有效位点计算
        if len(variant.ALT) > 1:
            gt[:] = GT_UNKNOWN
        chunk.append(gt)

        if len(chunk) >= chunk_size:
            yield np.array(chunk, dtype=np.int8)
            chunk = []

    if chunk:
        yield np.array(chunk, dtype=np.int8)

    vcf.close()


def compute_all_pairs_stats(
    gt_matrix: np.ndarray,
    pair_indices: np.ndarray,
) -> dict[str, np.ndarray]:
    """
    向量化计算所有样本对的统计量（无 Python 循环）。

    Args:
        gt_matrix: (variants, samples) int8 数组
        pair_indices: (n_pairs, 2) 样本索引对

    Returns:
        包含各统计量的字典，每个值为 (n_pairs,) 数组
    """
    # 提取所有样本对的基因型
    idx_a = pair_indices[:, 0]
    idx_b = pair_indices[:, 1]
    gt_a = gt_matrix[:, idx_a].T  # (n_pairs, n_variants)
    gt_b = gt_matrix[:, idx_b].T  # (n_pairs, n_variants)

    # 向量化计算（无 Python 循环）
    valid = (gt_a != GT_UNKNOWN) & (gt_b != GT_UNKNOWN)
    a_hom = (gt_a == GT_HOM_REF) | (gt_a == GT_HOM_ALT)
    b_hom = (gt_b == GT_HOM_REF) | (gt_b == GT_HOM_ALT)
    both_hom = a_hom & b_hom & valid

    # 基因型相等：gt_types 数值相同
    gt_equal = (gt_a == gt_b) & valid
    homo_equal = gt_equal & both_hom
    het_sites = ~both_hom & valid
    het_equal = gt_equal & het_sites

    # 沿 variants 轴求和
    return {
        "non_miss": valid.sum(axis=1),
        "a_hom": (a_hom & valid).sum(axis=1),
        "a_het": (~a_hom & valid).sum(axis=1),
        "b_hom": (b_hom & valid).sum(axis=1),
        "b_het": (~b_hom & valid).sum(axis=1),
        "both_hom": both_hom.sum(axis=1),
        "homo_equal": homo_equal.sum(axis=1),
        "het_sites": het_sites.sum(axis=1),
        "het_equal": het_equal.sum(axis=1),
    }


def format_results(
    pairs: list[tuple[str, str]],
    accumulators: np.ndarray,
) -> list[list[str | int | float]]:
    """
    将累加器数组转换为输出结果列表。

    Args:
        pairs: 样本对列表
        accumulators: (n_pairs, 10) 累加器数组
            列顺序: total_sites, non_miss, a_hom, a_het, b_hom, b_het,
                   both_hom, homo_equal, het_sites, het_equal

    Returns:
        结果列表，每个元素为一行输出
    """
    results = []
    for i, (sample_a, sample_b) in enumerate(pairs):
        total_sites = int(accumulators[i, 0])
        non_miss = int(accumulators[i, 1])
        a_hom = int(accumulators[i, 2])
        a_het = int(accumulators[i, 3])
        b_hom = int(accumulators[i, 4])
        b_het = int(accumulators[i, 5])
        both_hom = int(accumulators[i, 6])
        homo_equal = int(accumulators[i, 7])
        het_sites = int(accumulators[i, 8])
        het_equal = int(accumulators[i, 9])

        if non_miss == 0:
            results.append(
                [
                    sample_a,
                    sample_b,
                    total_sites,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0.0,
                    0,
                    0.0,
                    0,
                    0.0,
                    0,
                    0.0,
                ]
            )
            continue

        homo_pct = 0.0 if both_hom == 0 else round(100 * homo_equal / both_hom, 3)
        het_pct = 0.0 if het_sites == 0 else round(100 * het_equal / het_sites, 3)
        total_equal = homo_equal + het_equal
        total_equal_pct = round(100 * total_equal / non_miss, 3)
        total_diff = non_miss - total_equal
        total_diff_pct = round(100 * total_diff / non_miss, 3)

        results.append(
            [
                sample_a,
                sample_b,
                total_sites,
                non_miss,
                a_hom,
                a_het,
                b_hom,
                b_het,
                total_equal,
                total_equal_pct,
                homo_equal,
                homo_pct,
                het_equal,
                het_pct,
                total_diff,
                total_diff_pct,
            ]
        )

    return results


@app.command()
def main(
    vcf_file: Annotated[Path, typer.Argument(help="输入 VCF 文件")],
    output_file: Annotated[Path, typer.Argument(help="输出 CSV 文件")],
    compare_list: Annotated[
        Path | None,
        typer.Option("--compare-list", "-c", help="样本对比较列表文件（tab分隔，两列：A B）"),
    ] = None,
    chunk_size: Annotated[
        int,
        typer.Option("--chunk-size", "-s", help="分块大小（变异数/块），内存不足时可减小"),
    ] = 10000,
    verbose: Annotated[
        bool,
        typer.Option("--verbose", "-v", help="启用详细日志"),
    ] = False,
) -> None:
    """
    比较 VCF 文件中样本间的基因型一致性。

    计算所有样本对或指定样本对的相似度统计。
    输出包含纯合/杂合一致率等指标。
    采用 numpy 向量化计算，支持超大 VCF 文件和超多样本。
    """
    log_level = "DEBUG" if verbose else "INFO"
    logger.remove()
    logger.add(
        lambda msg: console.print(msg, end=""),
        level=log_level,
        format="<level>{message}</level>",
    )

    validate_file(vcf_file, "VCF 文件")
    if compare_list:
        validate_file(compare_list, "比较列表文件")

    # 获取样本列表
    samples = get_vcf_samples(vcf_file)
    sample_to_idx = {s: i for i, s in enumerate(samples)}
    logger.info(f"VCF 文件包含 {len(samples)} 个样本")

    # 构建比较对
    if compare_list:
        import polars as pl

        compare_df = pl.read_csv(compare_list, separator="\t", has_header=False, new_columns=["A", "B"])
        pairs = [(row[0], row[1]) for row in compare_df.iter_rows()]
        logger.info(f"从文件加载了 {len(pairs)} 个比较对")
    else:
        pairs = list(combinations(samples, 2))
        logger.info(f"将比较全部 {len(pairs)} 个样本对")

    # 过滤有效样本对并构建索引
    valid_pairs = [(a, b) for a, b in pairs if a in sample_to_idx and b in sample_to_idx]
    if len(valid_pairs) < len(pairs):
        skipped = len(pairs) - len(valid_pairs)
        logger.warning(f"跳过 {skipped} 个不存在的样本对")

    pair_indices = np.array([[sample_to_idx[a], sample_to_idx[b]] for a, b in valid_pairs], dtype=np.int32)

    logger.info(f"分块大小: {chunk_size}")

    # 初始化累加器数组
    # 列: total_sites, non_miss, a_hom, a_het, b_hom, b_het, both_hom, homo_equal, het_sites, het_equal
    accumulators = np.zeros((len(valid_pairs), 10), dtype=np.int64)

    # 分块处理
    chunk_count = 0
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("正在处理...", total=None)

        for gt_matrix in load_vcf_as_matrix(vcf_file, chunk_size):
            chunk_count += 1
            n_variants = gt_matrix.shape[0]
            progress.update(task, description=f"正在处理块 {chunk_count} ({n_variants} 变异)...")

            # 向量化计算当前块的统计量
            chunk_stats = compute_all_pairs_stats(gt_matrix, pair_indices)

            # 累加
            accumulators[:, 0] += n_variants  # total_sites
            accumulators[:, 1] += chunk_stats["non_miss"]
            accumulators[:, 2] += chunk_stats["a_hom"]
            accumulators[:, 3] += chunk_stats["a_het"]
            accumulators[:, 4] += chunk_stats["b_hom"]
            accumulators[:, 5] += chunk_stats["b_het"]
            accumulators[:, 6] += chunk_stats["both_hom"]
            accumulators[:, 7] += chunk_stats["homo_equal"]
            accumulators[:, 8] += chunk_stats["het_sites"]
            accumulators[:, 9] += chunk_stats["het_equal"]

    logger.info(f"共处理 {chunk_count} 个数据块，{int(accumulators[0, 0])} 个变异位点")

    # 格式化结果
    results = format_results(valid_pairs, accumulators)

    # 写入结果
    header = "A,B,总位点数,有效位点,A_纯合,A_杂合,B_纯合,B_杂合,整体相似度,整体相似度%,纯合相似度,纯合相似度%,杂合相似度,杂合相似度%,整体差异位点数,整体差异%\n"

    with open(output_file, "w", encoding="utf-8") as f:
        f.write(header)
        for row in results:
            line = ",".join(str(v) for v in row)
            f.write(f"{line}\n")

    console.print(f"[green]✓[/green] 结果已写入 {output_file}")


if __name__ == "__main__":
    app()
