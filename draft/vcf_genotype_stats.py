#!/usr/bin/env python3
"""
VCF基因型统计脚本 (修复版)

修复了缺失值 (./.) 被误判为纯合ALT的BUG，并优化了处理性能。
"""

from pathlib import Path
from typing import Dict, List
import typer
from cyvcf2 import VCF
from tqdm import tqdm


app = typer.Typer(add_completion=False, help="VCF基因型统计工具")


def get_genotype_class(a1: int, a2: int) -> str:
    """
    根据等位基因数值判断基因型分类 (直接处理整数，性能更高)

    Args:
        a1: 等位基因1 (cyvcf2中 -1 代表 .)
        a2: 等位基因2

    Returns:
        基因型分类: 'hom_ref', 'het', 'hom_alt', 'missing'
    """
    # cyvcf2 使用 -1 表示缺失值 (.)
    # 只要有一个等位基因缺失，通常视为该样本在该位点缺失或部分缺失
    if a1 == -1 or a2 == -1:
        return "missing"

    # 纯合REF (0/0)
    if a1 == 0 and a2 == 0:
        return "hom_ref"

    # 纯合ALT (1/1, 2/2 等)
    if a1 == a2:
        return "hom_alt"

    # 杂合型 (0/1, 1/2 等)
    return "het"


def process_vcf(vcf_path: Path, output_path: Path = None) -> Dict[str, int]:
    """
    处理VCF文件，统计基因型数量
    """
    vcf = VCF(str(vcf_path))

    # 统计总体基因型数量
    total_stats = {"hom_ref": 0, "het": 0, "hom_alt": 0, "missing": 0}

    # 结果缓存
    site_records = []

    print(f"处理VCF文件: {vcf_path}")
    print(f"检测到样本: {len(vcf.samples)}")
    print("开始处理位点...")

    # 遍历每个变异位点
    for variant in tqdm(vcf, desc="处理位点", unit="位点"):
        # 统计该位点所有样本的基因型
        site_hom_ref = 0
        site_het = 0
        site_hom_alt = 0
        site_missing = 0

        # variant.genotypes 返回列表: [[a1, a2, phased], ...]
        # 直接解包，不需要转字符串，速度快很多
        for gt in variant.genotypes:
            # 忽略倍性不是2的异常情况 (罕见)
            if len(gt) < 2:
                site_missing += 1
                continue

            gt_category = get_genotype_class(gt[0], gt[1])

            if gt_category == "hom_ref":
                site_hom_ref += 1
            elif gt_category == "het":
                site_het += 1
            elif gt_category == "hom_alt":
                site_hom_alt += 1
            else:
                site_missing += 1

        # 更新总体统计
        total_stats["hom_ref"] += site_hom_ref
        total_stats["het"] += site_het
        total_stats["hom_alt"] += site_hom_alt
        total_stats["missing"] += site_missing

        # 记录该位点的详细信息 (可选：只记录有变异的位点)
        if output_path:
            site_total = site_hom_ref + site_het + site_hom_alt + site_missing
            site_valid = site_hom_ref + site_het + site_hom_alt

            site_records.append(
                {
                    "CHROM": variant.CHROM,
                    "POS": variant.POS,
                    "REF": variant.REF,
                    "ALT": ",".join(variant.ALT) if variant.ALT else ".",
                    "hom_ref": site_hom_ref,
                    "het": site_het,
                    "hom_alt": site_hom_alt,
                    "missing": site_missing,
                    "total_samples": len(vcf.samples),
                    "het_rate": round(site_het / site_valid * 100, 2) if site_valid > 0 else 0,
                    "missing_rate": round(site_missing / site_total * 100, 2) if site_total > 0 else 0,
                }
            )

    vcf.close()

    # 输出结果到终端
    print("\n" + "=" * 60)
    print("统计结果 (修复版)")
    print("=" * 60)

    total_calls = sum(total_stats.values())
    valid_calls = total_stats["hom_ref"] + total_stats["het"] + total_stats["hom_alt"]

    print(f"\n总基因型调用数 (Samples * Sites): {total_calls:,}")

    print(f"\n分布情况:")
    print(f"  纯合REF (0/0): {total_stats['hom_ref']:,}")
    print(f"  杂合型 (0/1..): {total_stats['het']:,}")
    print(f"  纯合ALT (1/1..): {total_stats['hom_alt']:,}")
    print(f"  缺失/未检出 (./.): {total_stats['missing']:,}")

    if valid_calls > 0:
        het_rate = total_stats['het'] / valid_calls * 100
        missing_rate = total_stats['missing'] / total_calls * 100

        print(f"\n有效检出频率 (不含缺失):")
        print(f"  纯合REF: {total_stats['hom_ref']/valid_calls*100:.2f}%")
        print(f"  杂合型:  {total_stats['het']/valid_calls*100:.2f}%")
        print(f"  纯合ALT: {total_stats['hom_alt']/valid_calls*100:.2f}%")

        print(f"\n杂合率: {het_rate:.2f}%")
        print(f"缺失率: {missing_rate:.2f}%")

    # 输出到文件
    if output_path and site_records:
        import csv

        with open(output_path, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "CHROM",
                    "POS",
                    "REF",
                    "ALT",
                    "hom_ref",
                    "het",
                    "hom_alt",
                    "missing",
                    "total_samples",
                    "het_rate",
                    "missing_rate",
                ],
            )
            writer.writeheader()
            writer.writerows(site_records)
        print(f"\n详细结果已保存至: {output_path}")

    return total_stats


@app.command()
def analyze(
    vcf_file: Path = typer.Argument(..., help="输入VCF文件路径", exists=True),
    output_file: Path = typer.Option(None, "--output", "-o", help="输出CSV文件路径"),
):
    """
    分析VCF文件的基因型分布
    """
    try:
        process_vcf(vcf_file, output_file)
    except Exception as e:
        print(f"错误: {str(e)}")
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
