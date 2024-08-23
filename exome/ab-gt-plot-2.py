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


def gt2name(row, p1, p2, child, both_name):
    if row[child] == row[p1]:
        return p1
    if row[child] == row[p2]:
        return p2
    return both_name


def plot_origin(
    chr_df: pd.DataFrame,
    origin_df: pd.DataFrame,
    out_dir: Path,
    p1: str,
    p2: str,
    child_name: str,
    chr_file: Path,
    p1_color: str = "#E41A1C",
    p2_color: str = "#000080",
    both_color: str = "#696969",
    plot_both: bool = False,
    vertical: bool = False,
    both_name: str = "Both",
) -> None:
    origin_df = origin_df[origin_df[child_name] != "NN"].copy()
    origin_df["CHROM"] = origin_df["CHROM"].astype("str")
    origin_df = origin_df[origin_df["CHROM"].isin(chr_df.index.astype("str"))]

    my_gt2name = partial(gt2name, p1=p1, p2=p2, child=child_name, both_name=both_name)
    origin_df["origin"] = origin_df.apply(my_gt2name, axis=1)
    if not plot_both:
        origin_df = origin_df[origin_df["origin"] != both_name]

    color_map = {p1: p1_color, p2: p2_color, both_name: both_color}

    if vertical:
        out_table = out_dir / f"{child_name}.plot.tsv"
        out_plot_prefix = out_dir / f"{child_name}"
        origin_df["color"] = origin_df["origin"].map(color_map)
        origin_df.to_csv(out_table, sep="\t", index=False)
        plot_cmd = f'Rscript {VERTICAL_PLOT_R} -p {out_table} -c {chr_file} --p1_name {p1} --p2_name {p2} --p1_color "{p1_color}" --p2_color "{p2_color}" --both_name "{both_name}" --both_color "{both_color}" -o {out_plot_prefix}'
        print(plot_cmd)
        delegator.run(plot_cmd)
        out_table.unlink()
        return

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

    p1_patch = mpatches.Patch(color=p1_color, label=p1)
    p2_patch = mpatches.Patch(color=p2_color, label=p2)
    if plot_both:
        both_patch = mpatches.Patch(color=both_color, label=both_name)
        fig.legend(handles=[p1_patch, p2_patch, both_patch], ncol=3)
    else:
        fig.legend(handles=[p1_patch, p2_patch], ncol=2)
    out_png = out_dir / f"{child_name}.png"
    out_pdf = out_png.with_suffix(".pdf")
    plot_title = f"{child_name}"
    plt.title(plot_title)
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()


def main(
    gt_file: Path,
    chr_size: Path,
    out_dir: Path,
    p1: str,
    p2: str,
    child: Path,
    p1_color: str = "#0987ed",
    p2_color: str = "#0ec950",
    both_color: str = "#E41A1C",
    both_name: str = "Both",
    vertical: bool = False,
    plot_both: bool = False,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    child_list = [each.strip() for each in child.open()]
    chr_df = pd.read_table(
        chr_size, header=None, names=["chrom", "chrom_length"], index_col=0
    )
    gt_df = pd.read_table(gt_file)
    for child_name in tqdm(child_list):
        if child_name in [p1, p2]:
            child_df = gt_df[["CHROM", "POS", p1, p2]]
        else:
            child_df = gt_df[["CHROM", "POS", p1, p2, child_name]]
        plot_origin(
            chr_df,
            child_df,
            out_dir,
            p1,
            p2,
            child_name,
            chr_size,
            vertical=vertical,
            p1_color=p1_color,
            p2_color=p2_color,
            both_color=both_color,
            plot_both=plot_both,
            both_name=both_name,
        )


if __name__ == "__main__":
    typer.run(main)
