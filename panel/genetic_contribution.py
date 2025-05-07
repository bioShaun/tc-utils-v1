from functools import partial
from pathlib import Path
from typing import Optional

import delegator
import pandas as pd
import typer
from pysam import VariantFile

HOM_GT = [(0, 0), (1, 1)]

PLOT_R = Path(__file__).parent / "genetic_contribution.R"


def get_sample_indices(vcf_file, parent1, parent2, offspring):
    """获取样本在VCF文件中的索引"""
    vcf = VariantFile(vcf_file)
    samples = list(vcf.header.samples)

    try:
        p1_idx = samples.index(parent1)
        p2_idx = samples.index(parent2)
        off_idx = samples.index(offspring)
    except ValueError as e:
        print(f"错误: 样本名称不在VCF文件中: {e}")

    return p1_idx, p2_idx, off_idx


def filter_vcf(vcf_file: str, p1: str, p2: str, offspring: str) -> pd.DataFrame:
    p1_idx, p2_idx, off_idx = get_sample_indices(vcf_file, p1, p2, offspring)
    out_df_list = []
    for variant in VariantFile(vcf_file):
        alt = variant.alts

        if (alt is not None) and len(alt) > 1:
            continue
        gt_p1 = variant.samples[p1_idx]["GT"]
        gt_p2 = variant.samples[p2_idx]["GT"]
        gt_off = variant.samples[off_idx]["GT"]
        all_hom = all(gt in HOM_GT for gt in (gt_p1, gt_p2, gt_off))
        if not all_hom:
            continue
        out_df_list.append(
            {
                "chrom": variant.chrom,
                "pos": variant.pos,
                p1: gt_p1,
                p2: gt_p2,
                offspring: gt_off,
            }
        )
    return pd.DataFrame(out_df_list)


def cal_genetic_contribution(
    df_chr: pd.DataFrame, p1: str, p2: str, offspring: str, chrom: str
) -> dict:
    """计算遗传贡献率"""
    total_diff_count = df_chr.shape[0]
    p1_contribution = df_chr[df_chr[p1] == df_chr[offspring]].shape[0]
    p2_contribution = df_chr[df_chr[p2] == df_chr[offspring]].shape[0]
    p1_contribution_rate = 100 * p1_contribution / total_diff_count
    p2_contribution_rate = 100 * p2_contribution / total_diff_count
    return {
        "染色体": chrom,
        "总差异位点数": total_diff_count,
        f"{p1}差异位点数": p1_contribution,
        f"{p1}贡献率%": f"{p1_contribution_rate:.2f}",
        f"{p2}差异位点数": p2_contribution,
        f"{p2}贡献率%": f"{p2_contribution_rate:.2f}",
    }


def merge_genetic_contribution(
    df: pd.DataFrame, p1: str, p2: str, offspring: str, out_dir: Path
):
    """计算遗传贡献率"""
    filter_df = df[df[p1] != df[p2]]
    df_list = []
    for chrom, df_chr in filter_df.groupby("chrom"):
        df_list.append(cal_genetic_contribution(df_chr, p1, p2, offspring, chrom))
    if "chrom_group" in filter_df.columns:
        for chrom_group, df_chrom_group in filter_df.groupby("chrom_group"):
            df_list.append(
                cal_genetic_contribution(df_chrom_group, p1, p2, offspring, chrom_group)
            )
    df_list.append(cal_genetic_contribution(filter_df, p1, p2, offspring, "合计"))
    out_df = pd.DataFrame(df_list)
    out_df.to_excel(out_dir / f"{offspring}_genetic_contribution.xlsx", index=False)


def variant_site_class(row: pd.Series, p1: str, p2: str, offspring: str):
    """根据GT值判断变异类型"""
    if row[p1] == row[offspring] == row[p2]:
        return "all_same"
    if row[p1] == row[offspring]:
        return p1
    if row[p2] == row[offspring]:
        return p2
    return f"{offspring}_specific"


def plot(
    df: pd.DataFrame,
    p1: str,
    p2: str,
    offspring: str,
    chrom_info: Path,
    out_dir: Path,
    plot_bin=1.0,
    rscript_bin_path: Optional[Path] = None,
):
    """生成绘图数据"""
    plot_data_dir = out_dir / "plot_data"
    plot_data_dir.mkdir(parents=True, exist_ok=True)
    plot_df = df.copy()
    new_variant_site_class = partial(
        variant_site_class, p1=p1, p2=p2, offspring=offspring
    )
    plot_df["origin"] = plot_df.apply(new_variant_site_class, axis=1)
    plot_data = plot_data_dir / f"{offspring}_plot_data.tsv"
    plot_df.to_csv(plot_data, index=False, sep="\t")
    plot_dir = out_dir / "plot"
    plot_dir.mkdir(parents=True, exist_ok=True)
    if rscript_bin_path is not None:
        rscript_path = rscript_bin_path / "Rscript"
    else:
        rscript_path = "Rscript"
    plot_cmd = (
        f"{rscript_path} {PLOT_R} --plot_file {plot_data} --out_prefix {plot_dir}/{offspring} "
        f"--chr_size {chrom_info} --p1_name {p1} --p2_name {p2} --offspring_name {offspring}_specific --same_name all_same "
        f'--p1_color "yellow" --p2_color "#000080" --offspring_color "#E41A1C" --same_color "#696969" --bin_size {plot_bin}'
    )
    # print(plot_cmd)
    delegator.run(plot_cmd)


def main(
    vcf_file: str,
    parent_offspring_file: Path,
    chrom_info: Path,
    out_dir: Path,
    plot_bin: float = 1.0,
):
    out_dir.mkdir(parents=True, exist_ok=True)
    parent_offspring_df = pd.read_table(
        parent_offspring_file,
        header=None,
        names=[
            "offspring",
            "p1",
            "p2",
        ],
    )
    chrom_info_df = pd.read_table(chrom_info)
    for _, row in parent_offspring_df.iterrows():
        offspring, p1, p2 = row
        out_df = filter_vcf(vcf_file, p1, p2, offspring)
        out_df = out_df.merge(chrom_info_df)
        merge_genetic_contribution(out_df, p1, p2, offspring, out_dir)
        plot(out_df, p1, p2, offspring, chrom_info, out_dir, plot_bin=plot_bin)
        # out_df.to_csv(out_dir / f"{offspring}.csv", index=False)


if __name__ == "__main__":
    typer.run(main)
