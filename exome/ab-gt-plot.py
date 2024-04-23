from pathlib import Path
from typing import List, Tuple

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer
from tqdm import tqdm

COLOR_MAP = {"AA": "#E41A1C", "BB": "#000080", "AB": "#696969"}


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


def plot_origin(
    chr_df: pd.DataFrame,
    origin_df: pd.DataFrame,
    out_dir: Path,
    p1: str,
    p2: str,
    child_name: str,
) -> None:
    origin_df = origin_df.copy()
    origin_df.columns = ["Chr", "Pos", "Origin"]
    origin_df = origin_df[origin_df["Origin"] != "NN"]
    print(origin_df)
    plt.rcParams["font.size"] = 24

    # origin_df = pd.read_csv(origin_file)
    origin_df["Chr"] = origin_df["Chr"].astype("str")
    xaxis_ticks, xaxis_labels = get_plot_xaxis(chr_df)
    span = int(xaxis_ticks[-1] * 100000 / 50000000)
    fig, ax = plt.subplots()
    plot_width, plot_height = get_plot_size(chr_df)
    fig.set_figheight(plot_height)
    fig.set_figwidth(plot_width)

    for n, (chrom, chrom_length) in enumerate(chr_df.itertuples()):
        each_chrom_df = origin_df[origin_df["Chr"] == str(chrom)]
        y_pos = (plot_height - n - 2) * 10 + 4
        # ax.broken_barh([(0, chrom_length)], (y_pos, 6), facecolors=["#FFFFFF"])
        ax.broken_barh(
            [(0, chrom_length)],
            (y_pos, 6),
            facecolors=["#FFFFFF"],
            edgecolors=["#000000"],
        )
        if each_chrom_df.empty:
            continue

        plot_x = list(each_chrom_df.apply(lambda x: (int(x.Pos), span), axis=1))
        plot_colors = list(each_chrom_df["Origin"].map(lambda x: COLOR_MAP[x]))
        ax.broken_barh(plot_x, (y_pos, 6), facecolors=plot_colors)
    ax.set_yticks(
        [i * 10 + 7 for i in range(plot_height - 1)],
        labels=np.flip(chr_df.index),
    )
    ax.set_xticks(
        xaxis_ticks,
        xaxis_labels,
    )

    female_patch = mpatches.Patch(color="#E41A1C", label=p1)
    male_patch = mpatches.Patch(color="#000080", label=p2)
    both_patch = mpatches.Patch(color="#696969", label="BOTH")

    fig.legend(handles=[female_patch, male_patch, both_patch], ncol=3)
    out_png = out_dir / f"{child_name}.png"
    out_pdf = out_png.with_suffix(".pdf")
    plot_title = f"{child_name}"
    plt.title(plot_title)
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()


def main(
    gt_file: Path, chr_size: Path, out_dir: Path, p1: str, p2: str, child: Path
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    child_list = [each.strip() for each in child.open()]
    chr_df = pd.read_table(
        chr_size, header=None, names=["chrom", "chrom_length"], index_col=0
    )
    gt_df = pd.read_table(gt_file)
    for child_name in tqdm(child_list):
        child_df = gt_df[["CHROM", "POS", child_name]]
        plot_origin(chr_df, child_df, out_dir, p1, p2, child_name)


if __name__ == "__main__":
    typer.run(main)
