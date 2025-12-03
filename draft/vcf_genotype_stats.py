#!/usr/bin/env python3
"""
VCF基因型统计脚本

读取VCF文件，统计每个位点的基因型数量：
- 纯合REF型 (0/0)
- 杂合型 (0/1, 1/0, 0/2, 等)
- 纯合ALT型 (1/1, 2/2, 等)
"""

from pathlib import Path
from typing import Dict, Tuple
import typer
from cyvcf2 import VCF, Variant
from tqdm import tqdm


app = typer.Typer(add_completion=False, help="VCF基因型统计工具")


def parse_genotype(gt_string: str) -> str:
    """
    解析基因型字符串，返回分类结果

    Args:
        gt_string: 基因型字符串，如 '0/0', '0/1', '1/1', '1/2' 等

    Returns:
        基因型分类: 'hom_ref', 'het', 'hom_alt'
    """
    # 移除相位符号 (|, /)
    gt_parts = gt_string.replace('|', '/').split('/')

    if len(gt_parts) != 2:
        return 'unknown'

    try:
        a1, a2 = int(gt_parts[0]), int(gt_parts[1])
    except (ValueError, TypeError):
        return 'unknown'

    # 纯合REF (0/0)
    if a1 == 0 and a2 == 0:
        return 'hom_ref'

    # 纯合ALT (1/1, 2/2, 3/3, 等)
    elif a1 == a2 and a1 != 0:
        return 'hom_alt'

    # 杂合型 (0/1, 1/0, 0/2, 1/2, 等)
    else:
        return 'het'


def process_vcf(vcf_path: Path, output_path: Path = None) -> Dict[str, int]:
    """
    处理VCF文件，统计基因型数量

    Args:
        vcf_path: VCF文件路径
        output_path: 可选的输出文件路径

    Returns:
        包含统计结果的字典
    """
    vcf = VCF(str(vcf_path))

    # 统计总体基因型数量
    total_stats = {
        'hom_ref': 0,
        'het': 0,
        'hom_alt': 0,
        'unknown': 0
    }

    # 按染色体统计
    chrom_stats: Dict[str, Dict[str, int]] = {}

    # 按位置统计的详细记录
    site_records = []

    print(f"处理VCF文件: {vcf_path}")
    print(f"检测到样本: {len(vcf.samples)}")
    print("开始处理位点...")

    # 遍历每个变异位点
    for variant in tqdm(vcf, desc="处理位点", unit="位点"):
        chrom = variant.CHROM
        pos = variant.POS
        ref = variant.REF
        alt_alleles = variant.ALT  # ALT等位基因列表

        # 初始化染色体统计
        if chrom not in chrom_stats:
            chrom_stats[chrom] = {'hom_ref': 0, 'het': 0, 'hom_alt': 0, 'unknown': 0}

        # 统计该位点所有样本的基因型
        site_hom_ref = 0
        site_het = 0
        site_hom_alt = 0

        for sample_name in vcf.samples:
            gt = variant.genotypes[variant.sample_indexes[sample_name]]

            # gt是长度为3的列表: [allele1, allele2, phase]
            if len(gt) >= 2:
                gt_string = f"{gt[0]}/{gt[1]}"
            else:
                gt_string = str(gt[0]) if isinstance(gt[0], str) else "unknown"

            gt_category = parse_genotype(gt_string)

            # 统计该位点
            if gt_category == 'hom_ref':
                site_hom_ref += 1
            elif gt_category == 'het':
                site_het += 1
            elif gt_category == 'hom_alt':
                site_hom_alt += 1
            else:
                pass  # unknown类型不统计

            # 累加到总体统计
            total_stats[gt_category] += 1
            chrom_stats[chrom][gt_category] += 1

        # 记录该位点的详细信息
        if site_hom_ref + site_het + site_hom_alt > 0:  # 只记录有有效基因型的位点
            site_records.append({
                'CHROM': chrom,
                'POS': pos,
                'REF': ref,
                'ALT': ','.join(str(a) for a in alt_alleles) if alt_alleles else '.',
                'hom_ref': site_hom_ref,
                'het': site_het,
                'hom_alt': site_hom_alt,
                'total_samples': len(vcf.samples)
            })

    vcf.close()

    # 输出结果
    print("\n" + "="*60)
    print("统计结果")
    print("="*60)

    print(f"\n总位点数: {len(site_records)}")
    print(f"总样本数: {len(vcf.samples)}")

    print(f"\n总体基因型统计:")
    print(f"  纯合REF (0/0): {total_stats['hom_ref']:,}")
    print(f"  杂合型:         {total_stats['het']:,}")
    print(f"  纯合ALT:        {total_stats['hom_alt']:,}")
    if total_stats['unknown'] > 0:
        print(f"  未知/缺失:      {total_stats['unknown']:,}")

    total_called = total_stats['hom_ref'] + total_stats['het'] + total_stats['hom_alt']
    if total_called > 0:
        print(f"\n基因型频率:")
        print(f"  纯合REF: {total_stats['hom_ref']/total_called*100:.2f}%")
        print(f"  杂合型:  {total_stats['het']/total_called*100:.2f}%")
        print(f"  纯合ALT: {total_stats['hom_alt']/total_called*100:.2f}%")

    # 按染色体统计
    print(f"\n按染色体统计:")
    for chrom in sorted(chrom_stats.keys()):
        stats = chrom_stats[chrom]
        total = stats['hom_ref'] + stats['het'] + stats['hom_alt']
        if total > 0:
            print(f"  {chrom}: REF={stats['hom_ref']:,} | HET={stats['het']:,} | ALT={stats['hom_alt']:,}")

    # 输出到文件
    if output_path:
        import csv
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=[
                'CHROM', 'POS', 'REF', 'ALT', 'hom_ref', 'het', 'hom_alt', 'total_samples'
            ])
            writer.writeheader()
            writer.writerows(site_records)

        print(f"\n详细结果已保存至: {output_path}")

    return total_stats


@app.command()
def analyze(
    vcf_file: Path = typer.Argument(..., help="输入VCF文件路径", exists=True),
    output_file: Path = typer.Option(None, "--output", "-o", help="输出详细统计结果到CSV文件"),
    show_sites: bool = typer.Option(False, "--show-sites", help="显示前10个位点的详细信息"),
):
    """
    分析VCF文件的基因型分布
    """
    try:
        stats = process_vcf(vcf_file, output_file)

        if show_sites:
            print("\n前10个位点详情:")
            print("-" * 80)

    except FileNotFoundError:
        typer.echo(f"错误: 找不到文件 {vcf_file}", fg=typer.colors.RED)
        raise typer.Exit(1)
    except Exception as e:
        typer.echo(f"错误: 处理文件时发生异常 - {str(e)}", fg=typer.colors.RED)
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
