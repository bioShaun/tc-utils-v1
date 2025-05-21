#!/usr/bin/env python3
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
import pysam
import typer
from tqdm import tqdm


@dataclass
class AlleleStats:
    """存储等位基因统计信息的数据类"""

    missing_count: int = 0
    het_count: int = 0
    ref_count: int = 0
    alt_count: int = 0


@dataclass
class RateStats:
    """存储频率统计信息的数据类"""

    missing_rate: float = 0.0
    het_rate: float = 0.0
    maf: float = 0.0


def calculate_allele_stats(sample_genotypes: List[Any]) -> AlleleStats:
    """
    计算样本基因型的统计信息

    参数:
    sample_genotypes: 样本基因型列表

    返回:
    AlleleStats: 包含缺失样本数、杂合样本数、参考等位基因数和替代等位基因数的数据类
    """
    stats = AlleleStats()

    for sample in sample_genotypes:
        # 检查是否缺失
        if sample.get("GT") is None or all(g is None for g in sample["GT"]):
            stats.missing_count += 1
            continue

        # 获取非缺失的基因型
        genotype = [g for g in sample["GT"] if g is not None]

        # 检查是否杂合
        if len(set(genotype)) > 1:
            stats.het_count += 1

        # 统计参考和替代等位基因的数量
        for gt in genotype:
            if gt == 0:
                stats.ref_count += 1
            elif gt == 1:
                stats.alt_count += 1

    return stats


def calculate_rates(stats: AlleleStats, n_samples: int) -> RateStats:
    """
    计算缺失率、杂合率和MAF

    参数:
    stats: 等位基因统计信息
    n_samples: 总样本数

    返回:
    RateStats: 包含缺失率、杂合率和MAF的数据类
    """
    rate_stats = RateStats()

    # 计算缺失率
    rate_stats.missing_rate = stats.missing_count / n_samples if n_samples > 0 else 0

    # 计算杂合率
    non_missing_samples = n_samples - stats.missing_count
    rate_stats.het_rate = (
        stats.het_count / non_missing_samples if non_missing_samples > 0 else 0
    )

    # 计算MAF
    total_alleles = stats.ref_count + stats.alt_count
    if total_alleles > 0:
        ref_freq = stats.ref_count / total_alleles
        alt_freq = stats.alt_count / total_alleles
        rate_stats.maf = min(ref_freq, alt_freq)

    return rate_stats


def process_variant(record: Any) -> Optional[Dict]:
    """
    处理单个变异位点

    参数:
    record: VCF记录

    返回:
    Optional[Dict]: 包含变异位点统计信息的字典，如果不是双等位基因位点则返回None
    """
    # 跳过非双等位基因位点
    if len(record.alleles) != 2:
        return None

    # 获取基本信息
    chrom = record.chrom
    pos = record.pos
    alleles = ",".join(record.alleles)

    # 计算统计信息
    n_samples = len(record.samples)
    allele_stats = calculate_allele_stats(record.samples.values())
    rate_stats = calculate_rates(allele_stats, n_samples)

    # 返回结果
    return {
        "chrom": chrom,
        "pos": pos,
        "alleles": alleles,
        "missing": rate_stats.missing_rate,
        "het": rate_stats.het_rate,
        "maf": rate_stats.maf,
    }


def process_vcf(vcf_file: Path, output_file: Optional[Path] = None) -> pd.DataFrame:
    """
    从VCF文件中提取指定的信息并生成表格
    仅考虑2等位基因的情况（参考等位基因和一个替代等位基因）

    参数:
    vcf_file: VCF文件路径
    output_file: 输出表格文件路径，可选

    返回:
    DataFrame: 包含处理结果的DataFrame
    """
    # 打开VCF文件
    vcf = pysam.VariantFile(str(vcf_file))

    # 准备存储结果的列表
    results = []

    # 遍历VCF文件中的每个变异位点
    for record in tqdm(vcf, desc="Processing VCF"):
        result = process_variant(record)
        if result:
            results.append(result)

    # 创建DataFrame
    df = pd.DataFrame(results)

    # 保存到文件
    if output_file:
        df.to_csv(output_file, index=False)

    return df


def main(
    vcf_file: Path = typer.Argument(..., help="输入VCF文件路径"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="输出文件路径"),
):
    """
    从VCF文件提取信息并生成表格
    """
    # 如果未提供输出文件，则使用默认名称
    if output is None:
        output = Path(f"{vcf_file.stem}_stats.csv")

    # 处理VCF文件
    df = process_vcf(vcf_file, output)

    # 显示结果
    typer.echo(f"结果已保存到 {output}")
    typer.echo("\n数据预览:")
    typer.echo(df.head())


if __name__ == "__main__":
    typer.run(main)
