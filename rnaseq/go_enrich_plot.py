from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer
from matplotlib.patches import Rectangle


def plot_go(df: pd.DataFrame, compare: str, out_prefix: str) -> None:
    df["-log10(p.adjust)"] = -np.log10(df["p.adjust"])

    plot_width = 2 + int(0.5 * len(df))

    # 按类别分组并排序
    df_sorted = pd.concat(
        [
            df[df["Category"] == "biological_process"],
            df[df["Category"] == "cellular_component"],
            df[df["Category"] == "molecular_function"],
        ]
    )

    # 重置索引
    df_sorted = df_sorted.reset_index(drop=True)

    # 计算每个类别的数量
    bp_count = sum(df_sorted["Category"] == "biological_process")
    cc_count = sum(df_sorted["Category"] == "cellular_component")
    mf_count = sum(df_sorted["Category"] == "molecular_function")

    # 创建带间隔的x轴位置
    bar_width = 0.4  # 柱子宽度
    bar_spacing = 0.6  # 柱子之间的间距
    category_gap = 0.4  # 类别之间的间隔

    x_positions = []
    current_pos = 0

    # 为生物过程添加位置
    for i in range(bp_count):
        x_positions.append(current_pos)
        current_pos += bar_spacing
    current_pos += category_gap

    # 为细胞组分添加位置
    for i in range(cc_count):
        x_positions.append(current_pos)
        current_pos += bar_spacing
    current_pos += category_gap

    # 为分子功能添加位置
    for i in range(mf_count):
        x_positions.append(current_pos)
        current_pos += bar_spacing

    # 创建图形
    plt.figure(figsize=(plot_width, 6))

    # 设置颜色映射
    color_map = {
        "biological_process": "#E41A1C",  # 红色
        "cellular_component": "#377EB8",  # 蓝色
        "molecular_function": "#4DAF4A",  # 绿色
    }

    # 创建柱状图
    bars = plt.bar(x_positions, df_sorted["-log10(p.adjust)"], width=bar_width)

    # 为不同类别设置颜色
    for i, bar in enumerate(bars):
        bar.set_color(color_map[df_sorted["Category"].iloc[i]])

    # 设置坐标轴
    plt.xlabel("GO Terms", fontsize=12)
    plt.ylabel("-log10(p.adjust)", fontsize=12)
    plt.title(compare, fontsize=14)

    # 添加图例
    legend_elements = [
        Rectangle((0, 0), 1, 1, facecolor=color) for color in color_map.values()
    ]
    plt.legend(
        legend_elements,
        ["biological_process", "cellular_component", "molecular_function"],
        loc="upper right",
    )

    # 设置x轴标签
    # plt.xticks(x_positions, df_sorted["Description"], rotation=45, ha="right")
    plt.xticks(
        x_positions, [str(x) for x in df_sorted["Description"]], rotation=45, ha="right"
    )

    # 调整布局
    plt.tight_layout()

    plt.savefig(f"{out_prefix}.png", dpi=300)
    plt.savefig(f"{out_prefix}.pdf", format="pdf", bbox_inches="tight")


def main(enrich_file: Path, plot_count: int = 30) -> None:
    enrich_df = pd.read_csv(enrich_file)
    compare = enrich_file.stem.replace("_go_enrichment", "")
    plot_go(enrich_df[:plot_count], compare, enrich_file.stem)


if __name__ == "__main__":
    typer.run(main)
