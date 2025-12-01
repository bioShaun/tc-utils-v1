#!/usr/bin/env python3
"""
Select the highest-priority snpEff impact/effect per site and export to TSV.

Priority: HIGH > MODERATE > LOW > MODIFIER.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional

from cyvcf2 import VCF
import typer
from loguru import logger
from tqdm import tqdm


IMPACT_PRIORITY = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
IMPACT_RANK = {impact: rank for rank, impact in enumerate(IMPACT_PRIORITY)}


@dataclass
class AnnEntry:
    allele: str
    effect: str
    impact: str
    gene_id: str
    raw: str


def normalize_ann_values(raw_ann) -> List[str]:
    """Return ANN values as a flat list of strings."""
    if raw_ann is None:
        return []
    if isinstance(raw_ann, (list, tuple)):
        values = [str(x) for x in raw_ann]
    else:
        values = [str(raw_ann)]
    normalized: List[str] = []
    for item in values:
        normalized.extend([part for part in item.split(",") if part])
    return normalized


def parse_ann(entry: str) -> Optional[AnnEntry]:
    """Parse a single ANN entry into fields; return None if malformed."""
    parts = entry.split("|")
    if len(parts) < 3:
        return None
    allele, effect, impact = parts[0], parts[1], parts[2]
    gene_id = parts[4] if len(parts) > 4 else "."
    return AnnEntry(
        allele=allele or ".",
        effect=effect or ".",
        impact=impact or ".",
        gene_id=gene_id or ".",
        raw=entry,
    )


def pick_best_annotation(entries: Iterable[AnnEntry]) -> Optional[AnnEntry]:
    best: Optional[AnnEntry] = None
    best_rank = len(IMPACT_PRIORITY) + 1
    for entry in entries:
        rank = IMPACT_RANK.get(entry.impact.upper(), len(IMPACT_PRIORITY))
        if rank < best_rank:
            best = entry
            best_rank = rank
    return best


def extract_best_annotations(
    vcf_path: Path, output: Path, ann_tag: str = "ANN"
) -> None:
    vcf = VCF(str(vcf_path))
    output.parent.mkdir(parents=True, exist_ok=True)

    total = 0
    missing = 0

    with output.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(
            [
                "chrom",
                "pos",
                "effect",
                "impact",
                "gene_id",
            ]
        )

        for record in tqdm(vcf, desc="Extracting ANN", unit="record"):
            total += 1
            ann_values = normalize_ann_values(record.INFO.get(ann_tag))
            parsed = [ann for ann in (parse_ann(val) for val in ann_values) if ann]
            best = pick_best_annotation(parsed)
            if not best:
                missing += 1
                continue
            writer.writerow(
                [
                    record.CHROM,
                    record.POS,
                    best.effect,
                    best.impact,
                    best.gene_id,
                ]
            )

    logger.info(
        "Finished. Records processed: {}. Records without ANN: {}.",
        total,
        missing,
    )


def main(
    vcf: Path = typer.Argument(..., help="VCF/BCF file annotated by snpEff"),
    output: Path = typer.Option(
        Path("snpeff_impact.tsv"),
        "--output",
        "-o",
        help="Output TSV path",
    ),
    ann_tag: str = typer.Option(
        "ANN", "--ann-tag", "-t", help="INFO tag holding snpEff annotations"
    ),
):
    extract_best_annotations(vcf, output, ann_tag=ann_tag)


if __name__ == "__main__":
    typer.run(main)
