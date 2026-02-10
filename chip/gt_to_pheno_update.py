from pathlib import Path

import numpy as np
import pandas as pd
import typer

DISCLAIMER = (
    "本表基于当前已知的核心功能位点及基因单倍型信息，对相关性状进行注释。"
    "结果反映材料的遗传潜能，仅供分子辅助育种参考，实际田间表现受环境等多因素影响。"
)


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


def compute_gt_consistent(gt_allele_df: pd.DataFrame) -> pd.Series:
    ref = gt_allele_df["REF"].astype(str)
    alt = gt_allele_df["ALT"].astype(str)
    pheno_allele1 = gt_allele_df["ref_allele1"].astype(str)
    pheno_allele2 = gt_allele_df["ref_allele2"].astype(str)

    pheno_has_del = pheno_allele1.str.contains("del", na=False) | pheno_allele2.str.contains(
        "del", na=False
    )
    pheno_has_ins = pheno_allele1.str.contains("ins", na=False) | pheno_allele2.str.contains(
        "ins", na=False
    )

    ref_has_del = ref.str.contains("del", na=False)
    ref_has_ins = ref.str.contains("ins", na=False)
    alt_has_del = alt.str.contains("del", na=False)
    alt_has_ins = alt.str.contains("ins", na=False)

    ref_ok = ref.eq(".") | ref.eq(pheno_allele1) | ref.eq(pheno_allele2)
    alt_ok = alt.eq(".") | alt.eq(pheno_allele1) | alt.eq(pheno_allele2)
    result = ref_ok & alt_ok

    # Keep the original short-circuit order in filter_by_allele:
    # REF(del/ins) -> ALT(del/ins) -> fallback match check
    ref_del_mask = ref_has_del
    result = result.where(~ref_del_mask, pheno_has_del)

    ref_ins_mask = ~ref_del_mask & ref_has_ins
    result = result.where(~ref_ins_mask, pheno_has_ins)

    alt_del_mask = ~ref_del_mask & ~ref_ins_mask & alt_has_del
    result = result.where(~alt_del_mask, pheno_has_del)

    alt_ins_mask = ~ref_del_mask & ~ref_ins_mask & ~alt_del_mask & alt_has_ins
    result = result.where(~alt_ins_mask, pheno_has_ins)
    return result


def build_sample_phenotype_df(
    consistant_df: pd.DataFrame, sample_list: list[str]
) -> pd.DataFrame:
    sample_values = consistant_df[sample_list].to_numpy(dtype=object, copy=False)
    ref_allele1 = consistant_df["ref_allele1"].to_numpy()
    ref_allele2 = consistant_df["ref_allele2"].to_numpy()
    phenotype1 = consistant_df["phenotype1"].to_numpy()
    phenotype2 = consistant_df["phenotype2"].to_numpy()

    out = np.empty(sample_values.shape, dtype=object)
    for row_idx in range(sample_values.shape[0]):
        alleles = (ref_allele1[row_idx], ref_allele2[row_idx])
        phenos = (phenotype1[row_idx], phenotype2[row_idx])
        row_gts = sample_values[row_idx]
        unique_gts = pd.unique(row_gts)
        gt_map = {gt: get_gt_pheno(alleles, phenos, gt) for gt in unique_gts}
        out[row_idx] = [gt_map[gt] for gt in row_gts]
    return pd.DataFrame(out, columns=sample_list, index=consistant_df.index)


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
    gt_allele_df["gt_consistent"] = compute_gt_consistent(gt_allele_df)
    consistant_df = gt_allele_df[gt_allele_df["gt_consistent"]].copy()
    sample_list = gt_df.columns[4:].to_list()
    pheno_df = consistant_df[
        ["location_status", "CHROM", "POS", "category", "trait", "gene"]
    ].copy()
    sample_pheno_df = build_sample_phenotype_df(consistant_df, sample_list)
    pheno_df = pd.concat([pheno_df, sample_pheno_df], axis=1)
    out_df = format_output(pheno_df)
    with pd.ExcelWriter(out_file) as writer:
        pd.DataFrame([[DISCLAIMER]]).to_excel(
            writer, sheet_name="Sheet1", index=False, header=False
        )
        out_df.to_excel(writer, sheet_name="Sheet1", startrow=1)


if __name__ == "__main__":
    typer.run(main)
