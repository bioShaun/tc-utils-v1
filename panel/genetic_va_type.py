import re
from pathlib import Path

import delegator
import numpy as np
import pandas as pd
import typer
from loguru import logger

EXTRACT_VCF_FILEDS = (
    'CHROM POS REF ANN[*].ALLELE ANN[*].IMPACT "ANN[*].EFFECT" "ANN[*].HGVS_P"'
)

COLUMN_MAP = {
    "CHROM": "chrom",
    "POS": "pos",
    "REF": "refer",
    "ANN[*].ALLELE": "alt",
    "ANN[*].EFFECT": "effect",
    "ANN[*].IMPACT": "impact",
    "ANN[*].HGVS_P": "hgvs_p",
}

IMPACT_SCORE = {
    "HIGH": 0,
    "MODERATE": 1,
    "LOW": 2,
    "MODIFIER": 3,
}


def extract_snpeff_anno(vcf_path: Path, snpeff_dir: Path, force: bool) -> Path:
    annotation_file = vcf_path.with_suffix(".table.tsv.gz")
    if force or not annotation_file.is_file():
        extract_cmd = (
            f"zcat {vcf_path} | "
            f"{snpeff_dir}/scripts/vcfEffOnePerLine.pl | "
            f"java -jar {snpeff_dir}/SnpSift.jar extractFields - {EXTRACT_VCF_FILEDS} | "
            f"gzip > {annotation_file}"
        )
        logger.info(f"run: {extract_cmd}")
        delegator.run(extract_cmd)
    return annotation_file


def is_missense_or_nonsynonymous(hgvs_p: str) -> str:
    if pd.isna(hgvs_p):
        return "other"
    if "*" in hgvs_p or "?" in hgvs_p:
        return "non-synonymous"
    proteins = re.split("[0-9]+", hgvs_p.split(".")[1])
    if len(proteins) != 2:
        return "non-synonymous"
    protein_a, protein_b = proteins
    if protein_a == protein_b:
        return "synonymous"
    return "non-synonymous"


def main(
    snpeff_vcf: Path,
    out_file: Path,
    snpeff_dir: Path,
    force: bool = False,
) -> None:
    logger.info("extract snpeff annotation ...")
    snpeff_anno_file = extract_snpeff_anno(
        vcf_path=snpeff_vcf, snpeff_dir=snpeff_dir, force=force
    )
    snpeff_df = pd.read_table(snpeff_anno_file)
    snpeff_df.rename(columns=COLUMN_MAP, inplace=True)
    snpeff_df["alleles"] = snpeff_df["refer"].str.cat(snpeff_df["alt"], sep="/")
    snpeff_df["impact_score"] = snpeff_df["impact"].map(IMPACT_SCORE)
    snpeff_df.sort_values(by="impact_score", inplace=True)
    snpeff_df["effect"] = snpeff_df["effect"].map(lambda x: x.split("&")[0])
    snpeff_df.drop_duplicates(subset=["chrom", "pos", "alleles"], inplace=True)
    snpeff_df["variant_type"] = snpeff_df["hgvs_p"].map(is_missense_or_nonsynonymous)
    snpeff_df.drop(["impact_score", "hgvs_p", "refer", "alt"], axis=1, inplace=True)
    snpeff_df.to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)
