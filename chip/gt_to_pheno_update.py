from pathlib import Path

import pandas as pd
import typer


def filter_by_allele(pheno_alleles, gt_alleles):
    for gt_allele in gt_alleles:
        if "del" in gt_allele:
            return "del" in pheno_alleles
        if "ins" in gt_allele:
            return "ins" in pheno_alleles

    for gt_allele in gt_alleles:
        if gt_allele == ".":
            continue
        if gt_allele not in pheno_alleles:
            return False
    return True


def get_gt_pheno(alleles: tuple[str], phenos: tuple[str], gt: str) -> str:
    if "N" in gt:
        return "-"
    if "del" in alleles or "ins" in alleles:
        if "/" in gt:
            return "杂合"
        if "del" in alleles:
            return phenos[alleles.index("del")].capitalize()
        return phenos[alleles.index("ins")].capitalize()
    if len(set(gt)) > 1:
        return "杂合"
    return phenos[alleles.index(gt[0])].capitalize()


def gene_add_label(row):
    if row["location_status"] == "LC":
        return f'{row["gene"]}(LC)'
    return row["gene"]


def format_output(df: pd.DataFrame) -> pd.DataFrame:
    out_df = df.copy()
    out_df["gene"] = out_df.apply(gene_add_label, axis=1)
    rename_df = df.rename(columns={"trait": "性状", "gene": "基因"})
    rename_df = rename_df[rename_df.columns[4:]].set_index(["性状", "基因"])
    rename_df.columns.name = "样品"
    return rename_df.T


def main(gt_file: Path, pheno_file: Path, out_file: Path):
    gt_df = pd.read_table(gt_file)
    pheno_df = pd.read_table(pheno_file)
    gt_allele_df = pheno_df.rename(columns={"chrom": "CHROM", "pos": "POS"}).merge(
        gt_df
    )
    gt_allele_df["gt_consistent"] = gt_allele_df.apply(
        lambda row: filter_by_allele(
            (row["ref_allele1"], row["ref_allele2"]), (row["REF"], row["ALT"])
        ),
        axis=1,
    )
    consistant_df = gt_allele_df[gt_allele_df["gt_consistent"]].copy()
    sample_list = gt_df.columns[4:].to_list()
    pheno_df = consistant_df[
        ["location_status", "CHROM", "POS", "category", "trait", "gene"]
    ].copy()
    for sample in sample_list:
        pheno_df[sample] = consistant_df.apply(
            lambda row: get_gt_pheno(
                (row["ref_allele1"], row["ref_allele2"]),
                (row["phenotype1"], row["phenotype2"]),
                row[sample],
            ),
            axis=1,
        )
    out_df = format_output(pheno_df)
    out_df.to_excel(out_file, index=False)


if __name__ == "__main__":
    typer.run(main)
    typer.run(main)
