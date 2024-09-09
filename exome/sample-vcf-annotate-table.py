import csv
from pathlib import Path

import delegator
import pandas as pd
import typer

EXTRACT_VCF_FILEDS = (
    'CHROM POS REF "ANN[*].ALLELE" '
    '"ANN[*].IMPACT" "ANN[*].EFFECT" "ANN[*].GENEID"  '
    '"ANN[*].HGVS_C" "ANN[*].HGVS_P" GEN[0].GT GEN[0].AD'
)


COLUMN_MAP = {
    "ANN[*].ALLELE": "ALT",
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
        delegator.run(extract_cmd)


def main(
    vcf_file: Path,
    sample_name: str,
    snpeff_dir: Path,
    out_file: Path,
    force: bool = False,
):
    annotation_file = vcf_file.with_suffix(".raw.table.tsv.gz")
    extract_snpeff_anno(
        vcf_path=vcf_file,
        snpeff_dir=snpeff_dir,
        force=force,
        annotation_file=annotation_file,
    )
    df = pd.read_table(annotation_file)
    df.drop_duplicates(subset=["CHROM", "POS", "REF", "ALT"], inplace=True)
    df.rename(columns=COLUMN_MAP, inplace=True)
    df["Sample"] = sample_name
    df["Ref_Depth"] = df["AD"].map(lambda x: int(x[0]))
    df["Alt_Depth"] = df["AD"].map(lambda x: int(x[1]))
    if out_file.suffix != ".csv":
        raise ValueError(f"Output file should be .csv: {out_file}")
    df.to_csv(out_file, index=False, quoting=csv.QUOTE_NONNUMERIC)


if __name__ == "__main__":
    typer.run(main)
