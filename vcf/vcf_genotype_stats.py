#!/usr/bin/env python3
"""
VCF Genotype Statistics Tool

Calculates genotype frequencies (hom_ref, het, hom_alt, missing) across all samples in a VCF file.
Optimized for performance using cyvcf2.
"""

import csv
from pathlib import Path

import typer
from cyvcf2 import VCF
from loguru import logger
from tqdm import tqdm

app = typer.Typer(add_completion=False, help="VCF Genotype Statistics Tool")


def get_genotype_class(a1: int, a2: int) -> str:
    """
    Classify genotype based on allele values.

    Args:
        a1: Allele 1 value (-1 for missing, 0 for REF, >0 for ALT)
        a2: Allele 2 value

    Returns:
        One of: 'hom_ref', 'het', 'hom_alt', 'missing'
    """
    if a1 == -1 or a2 == -1:
        return "missing"
    if a1 == 0 and a2 == 0:
        return "hom_ref"
    if a1 == a2:
        return "hom_alt"
    return "het"


def process_vcf(vcf_path: Path, output_path: Path | None = None) -> dict[str, int]:
    """
    Process VCF file and aggregate genotype statistics.

    Args:
        vcf_path: Path to the input VCF file.
        output_path: Optional path to save site-level statistics as CSV.

    Returns:
        A dictionary containing aggregated genotype counts.
    """
    vcf = VCF(str(vcf_path))
    total_stats = {"hom_ref": 0, "het": 0, "hom_alt": 0, "missing": 0}
    site_records = []

    logger.info(f"Processing VCF: {vcf_path}")
    logger.info(f"Detected samples: {len(vcf.samples)}")

    for variant in tqdm(vcf, desc="Processing sites", unit="site"):
        site_stats = {"hom_ref": 0, "het": 0, "hom_alt": 0, "missing": 0}

        for gt in variant.genotypes:
            if len(gt) < 2:
                site_stats["missing"] += 1
                continue

            gt_cat = get_genotype_class(gt[0], gt[1])
            site_stats[gt_cat] += 1

        for key in total_stats:
            total_stats[key] += site_stats[key]

        if output_path:
            site_total = sum(site_stats.values())
            site_valid = site_total - site_stats["missing"]

            site_records.append({
                "CHROM": variant.CHROM,
                "POS": variant.POS,
                "REF": variant.REF,
                "ALT": ",".join(variant.ALT) if variant.ALT else ".",
                **site_stats,
                "total_samples": len(vcf.samples),
                "het_rate(%)": round(site_stats["het"] / site_valid * 100, 2) if site_valid > 0 else 0,
                "missing_rate(%)": round(site_stats["missing"] / site_total * 100, 2) if site_total > 0 else 0,
            })

    vcf.close()

    if output_path and site_records:
        save_site_records(output_path, site_records)

    return total_stats


def save_site_records(output_path: Path, records: list[dict]):
    """
    Save site-level genotype records to a CSV file.

    Args:
        output_path: Path to the output CSV file.
        records: List of dictionaries containing site statistics.
    """
    fieldnames = [
        "CHROM", "POS", "REF", "ALT", "hom_ref", "het", "hom_alt", "missing",
        "total_samples", "het_rate(%)", "missing_rate(%)"
    ]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
    logger.info(f"Detailed results saved to: {output_path}")


def print_summary(stats: dict[str, int]):
    """Print a summary of the genotype statistics."""
    total_calls = sum(stats.values())
    valid_calls = total_calls - stats["missing"]

    print("\n" + "=" * 60)
    print("Genotype Statistics Summary")
    print("=" * 60)
    print(f"\nTotal genotype calls (Samples * Sites): {total_calls:,}")
    print("\nDistribution:")
    print(f"  Homozygous REF (0/0):   {stats['hom_ref']:12,}")
    print(f"  Heterozygous (0/1..):   {stats['het']:12,}")
    print(f"  Homozygous ALT (1/1..): {stats['hom_alt']:12,}")
    print(f"  Missing (./.):          {stats['missing']:12,}")

    if valid_calls > 0:
        het_rate = stats['het'] / valid_calls * 100
        missing_rate = stats['missing'] / total_calls * 100

        print("\nFrequency (excluding missing):")
        print(f"  Homozygous REF: {stats['hom_ref']/valid_calls*100:6.2f}%")
        print(f"  Heterozygous:   {stats['het']/valid_calls*100:6.2f}%")
        print(f"  Homozygous ALT: {stats['hom_alt']/valid_calls*100:6.2f}%")
        print(f"\nHeterozygosity Rate: {het_rate:6.2f}%")
        print(f"Missingness Rate:    {missing_rate:6.2f}%")


@app.command()
def analyze(
    vcf_file: Path = typer.Argument(..., help="Path to input VCF file", exists=True, dir_okay=False),
    output: Path | None = typer.Option(None, "--output", "-o", help="Path to output CSV file for site-level stats"),
):
    """
    Analyze VCF genotype distribution.
    """
    try:
        stats = process_vcf(vcf_file, output)
        print_summary(stats)
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        raise typer.Exit(1) from e


if __name__ == "__main__":
    app()
