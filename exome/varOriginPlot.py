from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import typer


COLOR_MAP = {"Male_parent": "#377EB8", "Female_parent": "#4DAF4A", "Both": "#E41A1C"}


def get_plot_size(df: pd.DataFrame) -> Tuple[float, int]:
    height = len(df) + 1
    return 20, height


def get_plot_xaxis(df: pd.DataFrame) -> Tuple[List[int], List[str]]:
    chrom_max_length = df["chrom_length"].max()
    mega_base = np.floor(np.log10(chrom_max_length))
    max_chr_show_length = np.ceil(chrom_max_length / 10 ** mega_base)
    x_axis_ticks = [
        int(i * 10 ** mega_base) for i in range(int(max_chr_show_length) + 1)
    ]
    x_axis_labels = [
        f"{int(i * 10**mega_base / 1e6)}M" for i in range(int(max_chr_show_length) + 1)
    ]
    return x_axis_ticks, x_axis_labels


def main(chr_size: Path, origin_file: Path) -> None:
    plt.rcParams["font.size"] = 24
    chr_df = pd.read_table(
        chr_size, header=None, names=["chrom", "chrom_length"], index_col=0
    )
    origin_df = pd.read_csv(origin_file)
    origin_df["Chr"] = origin_df["Chr"].astype("str")
    xaxis_ticks, xaxis_labels = get_plot_xaxis(chr_df)
    span = xaxis_ticks[-1] * 100000 / 50000000
    fig, ax = plt.subplots()
    plot_width, plot_height = get_plot_size(chr_df)
    fig.set_figheight(plot_height)
    fig.set_figwidth(plot_width)

    for n, (chrom, chrom_length) in enumerate(chr_df.itertuples()):
        each_chrom_df = origin_df[origin_df["Chr"] == str(chrom)]
        plot_x = list(each_chrom_df.apply(lambda x: (x.Pos, span), axis=1))
        y_pos = (plot_height - n - 2) * 10 + 4
        plot_colors = list(each_chrom_df["Origin"].map(lambda x: COLOR_MAP[x]))
        ax.broken_barh([(0, chrom_length)], (y_pos, 6), facecolors=["#CFCFCF"])
        ax.broken_barh(plot_x, (y_pos, 6), facecolors=plot_colors)
    ax.set_yticks(
        [i * 10 + 7 for i in range(plot_height - 1)],
        labels=np.flip(chr_df.index),
    )
    ax.set_xticks(
        xaxis_ticks,
        xaxis_labels,
    )

    both_patch = mpatches.Patch(color="#E41A1C", label="Both")
    female_patch = mpatches.Patch(color="#377EB8", label="Male parent")
    male_patch = mpatches.Patch(color="#4DAF4A", label="Female parent")

    fig.legend(handles=[both_patch, female_patch, male_patch], ncol=3)
    plt.savefig(origin_file.with_suffix(".png"), dpi=300)
    plt.savefig(origin_file.with_suffix(".pdf"))


if __name__ == "__main__":
    typer.run(main)
