from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer

# 创建示例DataFrame


def create_dual_bar_chart(
    df: pd.DataFrame, figsize: Tuple[int, int] = (15, 8), dpi=300
):
    """
    创建双Y轴柱状图

    参数:
    df: DataFrame，包含 'Chromosome', 'Probe_Count', 'Length_Mb' 列
    figsize: 图表尺寸
    dpi: 分辨率
    """
    # 创建图形和主坐标轴
    fig, ax1 = plt.subplots(figsize=figsize, dpi=dpi)

    # 设置第一个Y轴（探针数量）
    x = np.arange(len(df))
    width = 0.35
    rects1 = ax1.bar(
        x - width / 2, df["Probe_Count"], width, label="Probe count", color="#4472C4"
    )
    ax1.set_ylabel("Probe count", fontsize=12)
    ax1.set_xlabel("Chromosome", fontsize=12)

    # 创建第二个Y轴（染色体长度）
    ax2 = ax1.twinx()
    rects2 = ax2.bar(
        x + width / 2,
        df["Length_Mb"],
        width,
        label="Chromosome length",
        color="#ED7D31",
    )
    ax2.set_ylabel("Chromosome length (Mb)", fontsize=12)

    # 设置X轴刻度
    ax1.set_xticks(x)
    ax1.set_xticklabels(df["Chromosome"])

    # 添加标题
    plt.title("Probe count vs. Chromosome length", fontsize=14, pad=20)

    # 添加网格线（只对第一个Y轴）
    ax1.grid(True, linestyle="--", alpha=0.3)

    # 合并两个坐标轴的图例
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right")

    # 调整布局
    plt.tight_layout()

    return fig


# 方法1：直接保存单张图表


def main(probe_bed: Path, chrom_size: Path, out_pdf: Path) -> None:
    chrom_len_df = pd.read_table(
        chrom_size, header=None, names=["Chromosome", "Length"]
    )
    chrom_len_df["Length_Mb"] = chrom_len_df["Length"] / 1e6
    probe_df = pd.read_table(
        probe_bed, header=None, names=["Chromosome", "Start", "End"]
    )

    df = pd.DataFrame(
        {
            "Chromosome": chrom_len_df["Chromosome"],
            "Probe_Count": probe_df.groupby("Chromosome").size(),
            "Length_Mb": chrom_len_df["Length_Mb"],
        }
    )
    figure_width = int(len(chrom_len_df) * 0.8)
    fig = create_dual_bar_chart(df, figsize=(figure_width, 8))
    fig.savefig(out_pdf, format="pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    typer.run(main)
