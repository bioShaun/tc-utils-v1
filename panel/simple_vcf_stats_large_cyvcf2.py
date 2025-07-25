#!/usr/bin/env python3
"""
VCF文件流式分析工具

此脚本使用cyvcf2库高效处理VCF文件，计算每个变异位点的统计信息。
支持流式处理大型文件，并提供并行计算功能。
"""

import csv
import multiprocessing
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import typer


@dataclass
class VcfStats:
    """存储VCF统计信息的数据类"""

    chrom: str
    pos: int
    alleles: str
    missing: float
    het: float
    maf: float

    def to_dict(self) -> Dict[str, Any]:
        """将数据类转换为字典"""
        return {
            "chrom": self.chrom,
            "pos": self.pos,
            "alleles": self.alleles,
            "missing": self.missing,
            "het": self.het,
            "maf": self.maf,
        }


@dataclass
class VariantData:
    """存储变异位点数据的可序列化数据类"""

    chrom: str
    pos: int
    ref: str
    alt: List[str]
    gt_types: List[int]
    aaf: float
    is_snp: bool


def variant_to_data(variant: Any) -> VariantData:
    """
    将cyvcf2.Variant对象转换为可序列化的VariantData

    参数:
    variant: cyvcf2.Variant对象

    返回:
    VariantData: 可序列化的变异位点数据
    """
    return VariantData(
        chrom=variant.CHROM,
        pos=variant.POS,
        ref=variant.REF,
        alt=list(variant.ALT),
        gt_types=variant.gt_types.tolist(),  # 将numpy数组转换为列表
        aaf=variant.aaf,
        is_snp=variant.is_snp,
    )


def process_variant(variant_data: VariantData) -> Optional[VcfStats]:
    """
    处理单个变异位点

    参数:
    variant_data: VariantData对象

    返回:
    Optional[VcfStats]: 如果是有效的SNP，返回统计信息，否则返回None
    """
    # 跳过非双等位基因位点
    if not variant_data.is_snp or len(variant_data.alt) != 1:
        return None

    # 获取基因型数据
    gt_types = np.array(variant_data.gt_types)  # 0=HOM_REF, 1=HET, 2=HOM_ALT, 3=UNKNOWN
    n_samples = len(gt_types)

    # 计算统计信息
    missing_count = np.sum(gt_types == 3)
    het_count = np.sum(gt_types == 1)

    # 计算等位基因频率
    alt_freq = variant_data.aaf
    ref_freq = 1 - alt_freq
    maf = min(ref_freq, alt_freq)

    # 计算缺失率和杂合率
    missing_rate = missing_count / n_samples
    non_missing = n_samples - missing_count
    het_rate = het_count / non_missing if non_missing > 0 else 0

    # 创建结果
    return VcfStats(
        chrom=variant_data.chrom,
        pos=variant_data.pos,
        alleles=f"{variant_data.ref},{','.join(variant_data.alt)}",
        missing=missing_rate,
        het=het_rate,
        maf=maf,
    )


def process_variant_batch(variant_data_list: List[VariantData]) -> List[VcfStats]:
    """
    处理一批变异位点

    参数:
    variant_data_list: 变异位点数据列表

    返回:
    List[VcfStats]: 有效SNP的统计信息列表
    """
    results = []
    for variant_data in variant_data_list:
        result = process_variant(variant_data)
        if result:
            results.append(result)
    return results


def write_results(writer: csv.DictWriter, results: List[VcfStats]) -> None:
    """
    将结果写入CSV文件

    参数:
    writer: CSV写入器
    results: 统计信息列表
    """
    for result in results:
        writer.writerow(result.to_dict())


def process_completed_futures(futures: Dict[Any, bool], writer: csv.DictWriter) -> int:
    """
    处理已完成的Future并写入结果

    参数:
    futures: Future字典
    writer: CSV写入器

    返回:
    int: 已处理的Future数量
    """
    done_futures = [f for f in list(futures.keys()) if f.done()]
    for future in done_futures:
        results = future.result()
        write_results(writer, results)
        futures.pop(future)
    return len(done_futures)


def open_vcf_file(vcf_file: Path):
    """
    打开VCF文件并返回cyvcf2.VCF对象

    参数:
    vcf_file: VCF文件路径

    返回:
    cyvcf2.VCF: VCF文件对象
    """
    try:
        import cyvcf2
    except ImportError:
        raise ImportError("请先安装cyvcf2库: pip install cyvcf2")

    return cyvcf2.VCF(str(vcf_file), gts012=True)


def print_progress(current: int, interval: int = 10000) -> None:
    """
    打印进度信息

    参数:
    current: 当前处理的记录数
    interval: 打印间隔
    """
    if current % interval == 0:
        sys.stdout.write(f"\r已处理 {current:,} 个变异位点...")
        sys.stdout.flush()


def stream_process_vcf(
    vcf_file: Path,
    output_file: Path,
    num_processes: Optional[int] = None,
    batch_size: int = 1000,
    show_progress: bool = True,
) -> Tuple[int, float]:
    """
    使用cyvcf2流式处理VCF文件

    参数:
    vcf_file: VCF文件路径
    output_file: 输出表格文件路径
    num_processes: 并行处理的进程数，默认为CPU核心数
    batch_size: 每批处理的记录数
    show_progress: 是否显示进度

    返回:
    Tuple[int, float]: (处理的变异位点数, 处理时间)
    """
    # 设置进程数
    if num_processes is None:
        num_processes = max(1, multiprocessing.cpu_count() - 1)

    start_time = time.time()

    # 打开VCF文件
    vcf = open_vcf_file(vcf_file)

    # 创建CSV写入器
    with open(output_file, "w", newline="") as f:
        fieldnames = ["chrom", "pos", "alleles", "missing", "het", "maf"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        # 使用ProcessPoolExecutor进行并行处理
        variant_count = process_variants_in_parallel(
            vcf, writer, num_processes, batch_size, show_progress
        )

    # 计算总处理时间
    elapsed_time = time.time() - start_time

    # 清理进度显示
    if show_progress:
        sys.stdout.write("\r" + " " * 50 + "\r")  # 清除当前行
        sys.stdout.flush()

    return variant_count, elapsed_time


def process_variants_in_parallel(
    vcf,
    writer: csv.DictWriter,
    num_processes: int,
    batch_size: int,
    show_progress: bool,
) -> int:
    """
    并行处理VCF文件中的变异位点

    参数:
    vcf: cyvcf2.VCF对象
    writer: CSV写入器
    num_processes: 并行处理的进程数
    batch_size: 每批处理的记录数
    show_progress: 是否显示进度

    返回:
    int: 处理的变异位点数
    """
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        # 初始化批次和futures
        batch = []
        futures = {}
        active_futures = 0
        max_active_futures = num_processes * 3  # 控制内存使用

        # 处理变异位点
        variant_count = 0
        for variant in vcf:
            # 转换为可序列化的数据
            variant_data = variant_to_data(variant)
            batch.append(variant_data)
            variant_count += 1

            if len(batch) >= batch_size:
                # 等待，如果活跃任务太多
                while active_futures >= max_active_futures:
                    processed = process_completed_futures(futures, writer)
                    if processed == 0:
                        time.sleep(0.1)  # 短暂等待
                    active_futures -= processed

                # 提交批次
                future = executor.submit(process_variant_batch, batch)
                futures[future] = True
                active_futures += 1
                batch = []

            # 更新进度
            if show_progress:
                print_progress(variant_count)

            # 定期检查完成的任务
            if variant_count % (batch_size * 10) == 0:
                active_futures -= process_completed_futures(futures, writer)

        # 处理最后一个批次
        if batch:
            future = executor.submit(process_variant_batch, batch)
            futures[future] = True

        # 等待所有任务完成
        while futures:
            active_futures -= process_completed_futures(futures, writer)
            if futures:  # 如果还有未完成的任务，短暂等待
                time.sleep(0.1)

    return variant_count


def main(
    vcf_file: Path = typer.Argument(..., help="输入VCF文件路径"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="输出文件路径"),
    processes: Optional[int] = typer.Option(
        None, "--processes", "-p", help="并行处理的进程数"
    ),
    batch_size: int = typer.Option(1000, "--batch-size", "-b", help="每批处理的记录数"),
    no_progress: bool = typer.Option(False, "--no-progress", help="不显示进度"),
):
    """
    从VCF文件提取信息并生成表格

    此脚本使用cyvcf2库进行高效的VCF文件处理，计算每个变异位点的缺失率、杂合率和次等位基因频率(MAF)。
    结果以TSV格式保存，可用于后续分析。

    示例:
        python vcf_analyzer.py input.vcf.gz -o output.tsv -p 4 -b 2000
    """
    # 如果未提供输出文件，则使用默认名称
    if output is None:
        output = Path(f"{vcf_file.stem}_stats.tsv")

    # 处理VCF文件
    try:
        typer.echo(f"开始处理VCF文件: {vcf_file}")
        typer.echo(f"使用进程数: {processes or '自动'}")
        typer.echo(f"批处理大小: {batch_size}")

        variant_count, elapsed_time = stream_process_vcf(
            vcf_file, output, processes, batch_size, not no_progress
        )

        # 显示结果
        typer.echo(f"处理完成!")
        typer.echo(f"总变异位点数: {variant_count}")
        typer.echo(f"处理时间: {elapsed_time:.2f}秒")
        typer.echo(f"处理速度: {variant_count/elapsed_time:.2f}变异位点/秒")
        typer.echo(f"结果已保存到: {output}")

        # 显示数据预览
        typer.echo("\n数据预览:")
        with open(output, "r") as f:
            for i, line in enumerate(f):
                if i == 0:
                    typer.echo(f"列名: {line.strip()}")
                elif i <= 5:
                    typer.echo(line.strip())
                else:
                    break

    except ImportError as e:
        typer.echo(f"错误: {e}", err=True)
        typer.echo("请安装cyvcf2: pip install cyvcf2", err=True)
        raise typer.Exit(code=1)
    except Exception as e:
        typer.echo(f"处理过程中发生错误: {e}", err=True)
        raise typer.Exit(code=1)


if __name__ == "__main__":
    typer.run(main)
