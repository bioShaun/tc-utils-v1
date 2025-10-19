import csv
from pathlib import Path
from typing import List, Tuple

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer
from cyvcf2 import VCF
from tqdm import tqdm

# ====================== 中文字体设置 ======================
# Ubuntu 常用中文字体: WenQuanYi Micro Hei / Noto Sans CJK SC / AR PL UKai CN
matplotlib.rcParams["font.sans-serif"] = ["AR PL UKai CN"]
matplotlib.rcParams["axes.unicode_minus"] = False  # 负号正常显示

# 颜色映射
COLOR_MAP = {
    "CONSISTENT": "#2ECC71",  # 绿色，与亲本一致
    "INCONSISTENT_HOM": "#7F7F7F",  # 灰色，不一致且纯和
    "INCONSISTENT_HET": "#E41A1C",  # 红色，不一致且杂合
}

VALID_CHILD_GT = {"0/0", "1/1", "0/1", "1/0"}  # 子代合法基因型


# ====================== VCF 提取 GT ======================
def read_vcf_gt(vcf_file: Path, samples: List[str]) -> pd.DataFrame:
    """
    使用 cyvcf2 读取 VCF 并提取 GT 信息，返回 DataFrame
    """
    vcf = VCF(str(vcf_file))
    sample_indices = [vcf.samples.index(s) for s in samples]

    rows = []
    for variant in vcf:
        chrom = variant.CHROM
        pos = variant.POS
        gts = []
        for idx in sample_indices:
            gt_tuple = variant.genotypes[idx]  # [allele1, allele2, phased, ...]
            if gt_tuple is None:
                gt_str = "./."
            else:
                gt_str = f"{gt_tuple[0]}/{gt_tuple[1]}"
            gts.append(gt_str)
        rows.append([chrom, pos] + gts)

    df = pd.DataFrame(rows, columns=["CHROM", "POS"] + samples)
    return df


# ====================== 绘图辅助 ======================
def get_plot_size(df: pd.DataFrame) -> Tuple[float, int]:
    height = len(df) + 1
    return 20, height


def get_plot_xaxis(df: pd.DataFrame) -> Tuple[List[int], List[str]]:
    chrom_max_length = df["chrom_length"].max()
    mega_base = np.floor(np.log10(chrom_max_length))
    max_chr_show_length = np.ceil(chrom_max_length / 10**mega_base)
    x_axis_ticks = [int(i * 10**mega_base) for i in range(int(max_chr_show_length) + 1)]
    x_axis_labels = [
        f"{int(i * 10**mega_base / 1e6)}M" for i in range(int(max_chr_show_length) + 1)
    ]
    return x_axis_ticks, x_axis_labels


def classify_origin(row: pd.Series, p1: str, child_name: str) -> str:
    """
    根据亲本与子代基因型判断类别
    """
    p1_gt = row[p1]
    child_gt = row[child_name]

    if child_gt in (".", "NN", "./."):
        return "NN"
    if child_gt == p1_gt:
        return "CONSISTENT"
    if child_gt in ("0/1", "1/0"):
        return "INCONSISTENT_HET"
    if child_gt in ("0/0", "1/1"):
        return "INCONSISTENT_HOM"
    return "NN"


# ====================== 作图 ======================
def plot_origin(
    chr_df: pd.DataFrame, gt_df: pd.DataFrame, out_dir: Path, p1: str, child_name: str
) -> dict:
    """
    作图并返回分类统计
    """
    origin_df = gt_df[["CHROM", "POS", p1, child_name]].copy()

    total_sites = len(origin_df)

    # 过滤亲本杂合
    origin_df = origin_df[origin_df[p1].isin(["0/0", "1/1"])]
    # 过滤子代非法基因型
    origin_df = origin_df[origin_df[child_name].isin(VALID_CHILD_GT)]

    # 分类
    origin_df["Origin"] = origin_df.apply(
        lambda x: classify_origin(x, p1, child_name), axis=1
    )
    origin_df = origin_df[origin_df["Origin"] != "NN"]

    # 分类统计
    counts = origin_df["Origin"].value_counts().to_dict()
    for key in ["CONSISTENT", "INCONSISTENT_HOM", "INCONSISTENT_HET"]:
        counts.setdefault(key, 0)
    total_retained = len(origin_df)
    counts["Total_retained"] = total_retained

    # 百分比
    counts["CONSISTENT_pct"] = (
        counts["CONSISTENT"] / total_retained * 100 if total_retained else 0
    )
    counts["INCONSISTENT_HOM_pct"] = (
        counts["INCONSISTENT_HOM"] / total_retained * 100 if total_retained else 0
    )
    counts["INCONSISTENT_HET_pct"] = (
        counts["INCONSISTENT_HET"] / total_retained * 100 if total_retained else 0
    )

    print(
        f"[{child_name}] 总位点: {total_sites}, 保留: {total_retained}, 绿色/灰色/红色计数: {counts}"
    )

    # ===== 作图 =====
    plt.rcParams["font.size"] = 24
    xaxis_ticks, xaxis_labels = get_plot_xaxis(chr_df)
    span = int(xaxis_ticks[-1] * 100000 / 50000000)
    fig, ax = plt.subplots()
    plot_width, plot_height = get_plot_size(chr_df)
    fig.set_figheight(plot_height)
    fig.set_figwidth(plot_width)

    for n, (chrom, chrom_length) in enumerate(chr_df.itertuples()):
        each_chrom_df = origin_df[origin_df["CHROM"] == str(chrom)]
        y_pos = (plot_height - n - 2) * 10 + 4
        ax.broken_barh(
            [(0, chrom_length)],
            (y_pos, 6),
            facecolors=["#FFFFFF"],
            edgecolors=["#000000"],
        )
        if each_chrom_df.empty:
            continue
        plot_x = list(each_chrom_df.apply(lambda x: (int(x.POS), span), axis=1))
        plot_colors = list(each_chrom_df["Origin"].map(lambda x: COLOR_MAP[x]))
        ax.broken_barh(plot_x, (y_pos, 6), facecolors=plot_colors)

    ax.set_yticks(
        [i * 10 + 7 for i in range(plot_height - 1)], labels=np.flip(chr_df.index)
    )
    ax.set_xticks(xaxis_ticks, xaxis_labels)

    legend_patches = [
        mpatches.Patch(color="#2ECC71", label=f"与 {p1} 一致"),
        mpatches.Patch(color="#7F7F7F", label="不一致纯和"),
        mpatches.Patch(color="#E41A1C", label="不一致杂合"),
    ]
    fig.legend(handles=legend_patches, ncol=3)

    out_png = out_dir / f"{child_name}.png"
    out_pdf = out_png.with_suffix(".pdf")
    plt.title(child_name)
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()

    return {"Child": child_name, **counts}


# ====================== 主函数 ======================
def main(vcf_file: Path, chr_size: Path, out_dir: Path, p1: str, child_file: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    child_list = [each.strip() for each in child_file.open()]
    chr_df = pd.read_table(
        chr_size, header=None, names=["chrom", "chrom_length"], index_col=0
    )

    # 使用 cyvcf2 提取 GT
    samples = [p1] + child_list
    gt_df = read_vcf_gt(vcf_file, samples)

    summary_list = []
    for child_name in tqdm(child_list):
        stats = plot_origin(chr_df, gt_df, out_dir, p1, child_name)
        summary_list.append(stats)

    # 保存汇总表
    summary_file = out_dir / "summary.csv"
    fieldnames = [
        "Child",
        "CONSISTENT",
        "INCONSISTENT_HOM",
        "INCONSISTENT_HET",
        "Total_retained",
        "CONSISTENT_pct",
        "INCONSISTENT_HOM_pct",
        "INCONSISTENT_HET_pct",
    ]
    with open(summary_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary_list:
            writer.writerow(row)

    print(f"汇总表已保存: {summary_file}")


# ====================== 入口 ======================
if __name__ == "__main__":
    typer.run(main)
