import csv
from pathlib import Path

import delegator
import pandas as pd
import typer
from loguru import logger

EXTRACT_VCF_FILEDS = (
    'CHROM POS REF "ALT" '
    '"ANN[*].IMPACT" "ANN[*].EFFECT" "ANN[*].GENEID"  '
    '"ANN[*].HGVS_C" "ANN[*].HGVS_P" GEN[0].GT GEN[0].AD'
)


COLUMN_MAP = {
    "ANN[*].IMPACT": "IMPACT",
    "ANN[*].EFFECT": "EFFECT",
    "ANN[*].GENEID": "GENEID",
    "ANN[*].HGVS_C": "HGVS_C",
    "ANN[*].HGVS_P": "HGVS_P",
    "GEN[0].GT": "GT",
    "GEN[0].AD": "AD",
}


def extract_snpeff_anno(
    vcf_path: Path, annotation_file: Path, snpeff_dir: Path, force: bool
) -> None:
    if force or not annotation_file.is_file():
        extract_cmd = (
            f"zcat {vcf_path} | "
            f"{snpeff_dir}/scripts/vcfEffOnePerLine.pl | "
            f"java -jar {snpeff_dir}/SnpSift.jar extractFields - {EXTRACT_VCF_FILEDS} | "
            f"gzip > {annotation_file}"
        )
        logger.info(f"run: {extract_cmd}")
        delegator.run(extract_cmd)


def get_sample_names(vcf_file: Path) -> list:
    cmd = f"bcftools query -l {vcf_file}"
    logger.info(f"run: {cmd}")
    return delegator.run(cmd).out.strip().split("\n")


def selectSampleVariants(
    vcf_file: Path, sample_name: str, out_file: Path, force: bool = False
) -> None:
    if force or not out_file.is_file():
        cmd = f"gatk SelectVariants -V {vcf_file} --exclude-non-variants --sample-name {sample_name} -O {out_file}"
        logger.info(f"run: {cmd}")
        delegator.run(cmd)


def main(
    vcf_file: Path,
    snpeff_dir: Path,
    out_dir: Path,
    force: bool = False,
):
    out_dir.mkdir(parents=True, exist_ok=True)
    sample_list = get_sample_names(vcf_file)
    for sample_name in sample_list:
        sample_vcf = vcf_file.with_suffix(f".{sample_name}.vcf.gz")
        selectSampleVariants(
            vcf_file=vcf_file,
            sample_name=sample_name,
            out_file=sample_vcf,
            force=force,
        )
        annotation_file = sample_vcf.with_suffix(".raw.table.tsv.gz")
        extract_snpeff_anno(
            vcf_path=sample_vcf,
            snpeff_dir=snpeff_dir,
            force=force,
            annotation_file=annotation_file,
        )
        df = pd.read_table(annotation_file)
        df.rename(columns=COLUMN_MAP, inplace=True)
        df.drop_duplicates(subset=["CHROM", "POS", "REF", "HGVS_C"], inplace=True)
        df["Sample"] = sample_name
        df["Ref_Depth"] = df["AD"].map(lambda x: int(x[0]))
        df["Alt_Depth"] = df["AD"].map(lambda x: ",".join(x[1]))
        sample_vcf_table = out_dir / f"{sample_name}.variant.table.csv"
        df.to_csv(sample_vcf_table, index=False, quoting=csv.QUOTE_NONNUMERIC)


if __name__ == "__main__":
    typer.run(main)
