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
        return "Missing"
    if "del" in alleles or "ins" in alleles:
        if "/" in gt:
            return "Heterozygous"
        if "del" in alleles:
            return phenos[alleles.index("del")].capitalize()
        return phenos[alleles.index("ins")].capitalize()
    if len(set(gt)) > 1:
        return "Heterozygous"
    return phenos[alleles.index(gt[0])].capitalize()


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
    sample_list = consistant_df.columns[10:-1].to_list()
    pheno_df = consistant_df[consistant_df.columns[:4]][
        ["CHROM", "POS", "trait", "gene"]
    ].copy()
    # 用字典先保存每个 sample 的新列数据
    sample_pheno_data = {}

    for sample in sample_list:
        sample_pheno_data[sample] = consistant_df.apply(
            lambda row: get_gt_pheno(
                (row["ref_allele1"], row["ref_allele2"]),
                (row["phenotype1"], row["phenotype2"]),
                row[sample],
            ),
            axis=1,
        )

    # 一次性拼接成新的 DataFrame
    pheno_sample_df = pd.DataFrame(sample_pheno_data)

    # 合并到原始 pheno_df，如果 pheno_df 是空的，可以直接用这个新 df
    pheno_df = pd.concat([pheno_df, pheno_sample_df], axis=1)
    pheno_df.to_excel(out_file, index=False)


if __name__ == "__main__":
    typer.run(main)
