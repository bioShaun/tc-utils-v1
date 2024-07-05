from functools import partial
from pathlib import Path
from typing import List, Tuple

import delegator
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer
from tqdm import tqdm

COLOR_MAP = {"AA": "#E41A1C", "BB": "#000080", "AB": "#696969"}
VERTICAL_PLOT_R = Path(__file__).parent / "ab-gt-plot-2-vertical-plot.R"


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


def gt2name(gt: str):
    if gt in ["./.", "."]:
        return "MISS"
    if gt == "0/0":
        return "REF"
    if "0" in gt:
        return "HET"
    return "ALT"


def plot_gt(
    chr_df: pd.DataFrame,
    origin_df: pd.DataFrame,
    sample_name: str,
    out_dir: Path,
    ref_color: str = "#E41A1C",
    alt_color: str = "#000080",
    het_color: str = "#696969",
    miss_color: str = "#FFFFFF",
) -> None:
    origin_df["CHROM"] = origin_df["CHROM"].astype("str")
    origin_df = origin_df[origin_df["CHROM"].isin(chr_df.index.astype("str"))]

    origin_df["genotype"] = origin_df[sample_name].map(gt2name)

    color_map = {
        "REF": ref_color,
        "HET": het_color,
        "ALT": alt_color,
        "MISS": miss_color,
    }

    plt.rcParams["font.size"] = 24

    # origin_df = pd.read_csv(origin_file)
    xaxis_ticks, xaxis_labels = get_plot_xaxis(chr_df)
    span = int(xaxis_ticks[-1] * 100000 / 50000000)
    fig, ax = plt.subplots()

    plot_width, plot_height = get_plot_size(chr_df)
    fig.set_figheight(plot_height)
    fig.set_figwidth(plot_width)

    for n, (chrom, chrom_length) in enumerate(chr_df.itertuples()):
        each_chrom_df = origin_df[origin_df["CHROM"] == str(chrom)]
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

        plot_x = list(each_chrom_df.apply(lambda x: (int(x.POS), span), axis=1))
        plot_colors = list(each_chrom_df["origin"].map(lambda x: color_map[x]))
        ax.broken_barh(plot_x, (y_pos, 6), facecolors=plot_colors)
    ax.set_yticks(
        [i * 10 + 7 for i in range(plot_height - 1)],
        labels=np.flip(chr_df.index),
    )
    ax.set_xticks(
        xaxis_ticks,
        xaxis_labels,
    )

    ref_patch = mpatches.Patch(color=ref_color, label="REF")
    alt_patch = mpatches.Patch(color=alt_color, label="ALT")
    het_patch = mpatches.Patch(color=het_color, label="HET")
    miss_patch = mpatches.Patch(color=miss_color, label="MISS")
    fig.legend(handles=[ref_patch, alt_patch, het_patch, miss_patch], ncol=4)
    out_png = out_dir / f"{sample_name}.png"
    out_pdf = out_png.with_suffix(".pdf")
    plot_title = f"{sample_name}"
    plt.title(plot_title)
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()


def main(
    gt_file: Path,
    chr_size: Path,
    out_dir: Path,
    ref_color: str = "#0987ed",
    alt_color: str = "#0ec950",
    het_color: str = "#E41A1C",
    miss_color: str = "#E41A1C",
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    chr_df = pd.read_table(
        chr_size, header=None, names=["chrom", "chrom_length"], index_col=0
    )
    gt_df = pd.read_table(gt_file)
    for sample_name in tqdm(gt_df.columns[4:]):
        sample_df = gt_df[["CHROM", "POS", sample_name]]
        plot_gt(
            chr_df,
            sample_df,
            sample_name,
            out_dir,
            ref_color=ref_color,
            alt_color=alt_color,
            het_color=het_color,
            miss_color=miss_color,
        )


if __name__ == "__main__":
    typer.run(main)
