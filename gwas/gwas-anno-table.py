from pathlib import Path

import delegator
import pandas as pd
import typer

EXTRACT_VCF_FILEDS = 'CHROM POS REF ALT "ANN[*].EFFECT" "ANN[*].FEATUREID" "ANN[*].GENEID"  "ANN[*].HGVS_C" "ANN[*].HGVS_P"'


COLUMN_MAP = {
    "CHROM": "chrom",
    "POS": "pos",
    "REF": "ref",
    "ALT": "alt",
    "ANN[*].EFFECT": "effect",
    "ANN[*].FEATUREID": "transcript_id",
    "ANN[*].GENEID": "gene_id",
    "ANN[*].HGVS_C": "hgvs_c",
    "ANN[*].HGVS_P": "hgvs_p",
}


def extract_snpeff_anno(vcf_path: Path, snpeff_path: Path) -> Path:
    annotation_file = vcf_path.with_suffix(".table.tsv.gz")
    extract_ann_cmd = (
        f"zcat {vcf_path} | {snpeff_path}/scripts/vcfEffOnePerLine.pl | "
        f"java -jar {snpeff_path}/SnpSift.jar extractFields - {EXTRACT_VCF_FILEDS} | gzip> {annotation_file}"
    )
    delegator.run(extract_ann_cmd)
    return annotation_file


def main(
    snpeff_vcf: Path, eggnog_anno: Path, out_file: Path, snpeff_path: Path
) -> None:
    snpeff_anno_file = extract_snpeff_anno(vcf_path=snpeff_vcf, snpeff_path=snpeff_path)
    eggnog_df = pd.read_table(eggnog_anno, sep="\t", skiprows=4, usecols=[0, 7])
    eggnog_df.columns = ["transcript_id", "description"]
    snpeff_df = pd.read_table(snpeff_anno_file)
    snpeff_df.drop_duplicates(subset=["chrom", "pos"])
    merged_df = snpeff_df.merge(eggnog_df, how="left")
    merged_df.fillna("--", inplace=True)
    merged_df.to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)
