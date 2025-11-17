#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
根据输入的位点ID列表，在VCF中寻找邻近的位点，并按MAF挑选替代位点。

输入:
    - 一个包含形如 chrom_pos 的ID列表文件（以空白分隔取首列）。
    - 需要搜索的VCF/BCF文件。

逻辑:
    1) 读取VCF，计算每个位点的MAF，并按染色体存入有序列表。
    2) 对每个目标位点，在同染色体上按距离选出最近的N个邻近位点（默认跳过同位置的位点）。
    3) 从这N个候选中，取MAF最高的m个作为替换并输出。
"""

from __future__ import annotations

import csv
from bisect import bisect_left
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from cyvcf2 import VCF
import typer
from loguru import logger
from tqdm import tqdm


@dataclass
class VariantSummary:
    chrom: str
    pos: int
    vid: str
    maf: float
    ref: str
    alt: str


@dataclass
class ChromVariants:
    variants: List[VariantSummary]
    positions: List[int]


@dataclass
class TargetSite:
    chrom: str
    pos: int
    raw_id: str


def parse_target_ids(id_file: Path) -> List[TargetSite]:
    targets: List[TargetSite] = []
    with id_file.open() as f:
        for line_no, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            site_id = line.split()[0]
            try:
                chrom, pos_str = site_id.split("_", 1)
                pos = int(pos_str)
            except ValueError as err:
                raise ValueError(
                    f"无法解析第 {line_no} 行的ID `{site_id}`，应为 chrom_pos 格式"
                ) from err
            targets.append(TargetSite(chrom=chrom, pos=pos, raw_id=site_id))
    if not targets:
        raise ValueError("ID列表为空，无法继续")
    return targets


def compute_genotype_metrics_simple(record) -> Dict[str, Optional[float]]:
    """使用 cyvcf2 内置属性的简化版本计算 MAF/缺失/杂合率。"""
    total_samples = len(record.genotypes)

    # 使用内置属性
    called_samples = record.num_het + record.num_hom_ref + record.num_hom_alt
    missing_samples = record.num_unknown
    het_samples = record.num_het

    # 计算等位基因计数
    num_alleles = len(record.ALT) + 1
    allele_counts = [0] * num_alleles

    for gt in record.genotypes:
        a1, a2 = gt[0], gt[1]
        if a1 >= 0:
            allele_counts[a1] += 1
        if a2 >= 0:
            allele_counts[a2] += 1

    total_called_alleles = sum(allele_counts)

    # MAF
    maf: Optional[float] = None
    if total_called_alleles > 0:
        freqs = [c / total_called_alleles for c in allele_counts if c > 0]
        if len(freqs) > 1:
            maf = min(freqs)

    missing_rate = missing_samples / total_samples if total_samples > 0 else None
    het_rate = het_samples / called_samples if called_samples > 0 else None

    return {
        "maf": maf,
        "missing_rate": missing_rate,
        "het_rate": het_rate,
        "called_samples": called_samples,
        "total_samples": total_samples,
    }


def calculate_maf(record) -> float:
    metrics = compute_genotype_metrics_simple(record)
    maf = metrics["maf"]
    if maf is not None:
        return maf

    # 如果基因型缺失，尝试使用 INFO/AF 回退
    info_af = record.INFO.get("AF")
    alt_freqs: Sequence[float] = info_af if info_af else ()
    if alt_freqs:
        ref_freq = max(0.0, 1.0 - sum(alt_freqs))
        freqs = [ref_freq, *alt_freqs]
        return min(freqs)
    return 0.0


def build_variant_index(vcf_path: Path) -> Dict[str, ChromVariants]:
    vcf = VCF(str(vcf_path))
    variants_by_chrom: Dict[str, List[VariantSummary]] = {}

    for record in tqdm(vcf, desc="读取VCF并计算MAF"):
        maf = calculate_maf(record)
        vid = record.ID or f"{record.CHROM}_{record.POS}"
        alt = ",".join(record.ALT) if record.ALT else ""
        variants_by_chrom.setdefault(record.CHROM, []).append(
            VariantSummary(
                chrom=record.CHROM,
                pos=record.POS,
                vid=vid,
                maf=maf,
                ref=record.REF,
                alt=alt,
            )
        )

    index: Dict[str, ChromVariants] = {}
    for chrom, variants in variants_by_chrom.items():
        variants.sort(key=lambda v: v.pos)
        positions = [v.pos for v in variants]
        index[chrom] = ChromVariants(variants=variants, positions=positions)

    if not index:
        raise ValueError("未能从VCF读取到任何位点")
    return index


def nearest_neighbors(
    chrom_variants: ChromVariants, pos: int, neighbor_count: int
) -> List[VariantSummary]:
    positions = chrom_variants.positions
    variants = chrom_variants.variants
    idx = bisect_left(positions, pos)
    left = idx - 1
    right = idx
    neighbors: List[VariantSummary] = []

    while len(neighbors) < neighbor_count and (left >= 0 or right < len(variants)):
        left_dist = pos - positions[left] if left >= 0 else None
        right_dist = positions[right] - pos if right < len(variants) else None

        if right_dist is None or (left_dist is not None and left_dist <= right_dist):
            neighbors.append(variants[left])
            left -= 1
        else:
            neighbors.append(variants[right])
            right += 1

    return neighbors


def select_replacements(
    candidates: Iterable[VariantSummary], pos: int, top_m: int, allow_same_pos: bool
) -> List[VariantSummary]:
    filtered = []
    for c in candidates:
        if not allow_same_pos and c.pos == pos:
            continue
        filtered.append(c)

    sorted_candidates = sorted(
        filtered, key=lambda v: (-v.maf, abs(v.pos - pos), v.pos)
    )
    return sorted_candidates[:top_m]


def main(
    id_list: Path = typer.Argument(..., help="包含 chrom_pos ID 的文本文件"),
    vcf: Path = typer.Argument(..., help="输入VCF/BCF文件"),
    output: Path = typer.Option(
        Path("vcf_site_replace.tsv"), "--output", "-o", help="输出TSV路径"
    ),
    neighbors: int = typer.Option(
        10,
        "--neighbors",
        "-n",
        help="按距离最近的顺序挑选的候选数（同染色体前后位点合计N个）",
    ),
    top: int = typer.Option(3, "--top", "-m", help="从候选中选出的位点数"),
    allow_same_pos: bool = typer.Option(
        False, help="是否允许选择与目标位置相同的位点作为替换"
    ),
):
    if neighbors <= 0:
        raise typer.BadParameter("--neighbors 必须大于0")
    if top <= 0:
        raise typer.BadParameter("--top 必须大于0")
    if top > neighbors:
        logger.warning("--top 大于 --neighbors，实际只会返回可用的候选数量")

    targets = parse_target_ids(id_list)
    index = build_variant_index(vcf)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(
            [
                "target_id",
                "target_chrom",
                "target_pos",
                "replacement_rank",
                "replacement_id",
                "replacement_chrom",
                "replacement_pos",
                "maf",
                "distance",
                "ref",
                "alt",
            ]
        )

        for target in targets:
            chrom_data = index.get(target.chrom)
            if not chrom_data:
                logger.warning(f"{target.raw_id} 未在VCF中找到同染色体数据，跳过")
                continue

            neighbors_for_target = nearest_neighbors(
                chrom_data, target.pos, neighbor_count=neighbors
            )
            replacements = select_replacements(
                neighbors_for_target,
                pos=target.pos,
                top_m=top,
                allow_same_pos=allow_same_pos,
            )
            if not replacements:
                logger.warning(f"{target.raw_id} 附近未找到可替换位点")
                continue

            for rank, variant in enumerate(replacements, start=1):
                writer.writerow(
                    [
                        target.raw_id,
                        target.chrom,
                        target.pos,
                        rank,
                        variant.vid,
                        variant.chrom,
                        variant.pos,
                        f"{variant.maf:.5f}",
                        abs(variant.pos - target.pos),
                        variant.ref,
                        variant.alt,
                    ]
                )

    logger.info(f"完成，结果已写入 {output}")


if __name__ == "__main__":
    typer.run(main)
