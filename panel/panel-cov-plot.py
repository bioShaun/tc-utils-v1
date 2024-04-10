from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer
from loguru import logger
from tqdm import tqdm

COLOR_MAP = {"AA": "#E41A1C", "BB": "#000080", "AB": "#696969"}


def merge_chr(df: pd.DataFrame, split_bed: Path) -> pd.DataFrame:
    split_bed_df = pd.read_csv(
        split_bed,
        header=None,
        names=["new_chrom", "offset", "offset_end", "chrom"],
        sep="\t",
    )
    print(df)
    print(split_bed_df)
    merged_df = df.merge(split_bed_df)
    merged_df["new_start"] = merged_df["start"] + merged_df["offset"]
    merged_df["new_end"] = merged_df["end"] + merged_df["offset"]
    merged_df.drop(
        ["chrom", "offset", "offset_end", "start", "end"], axis=1, inplace=True
    )
    merged_df.rename(
        columns={"new_chrom": "chrom", "new_start": "start", "new_end": "end"},
        inplace=True,
    )
    return merged_df


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


def plot_cov(
    cov_df: pd.DataFrame, chr_df: pd.DataFrame, out_dir: Path, name: str
) -> None:
    plt.rcParams["font.size"] = 24
    cov_df["chrom"] = cov_df["chrom"].astype("str")
    xaxis_ticks, xaxis_labels = get_plot_xaxis(chr_df)
    span = int(xaxis_ticks[-1] * 100000 / 50000000)
    fig, ax = plt.subplots()
    plot_width, plot_height = get_plot_size(chr_df)
    fig.set_figheight(plot_height)
    fig.set_figwidth(plot_width)

    for n, (chrom, chrom_length) in enumerate(chr_df.itertuples()):
        each_chrom_df = cov_df[cov_df["chrom"] == str(chrom)]
        y_pos = (plot_height - n - 2) * 10 + 4
        ax.broken_barh(
            [(0, chrom_length)],
            (y_pos, 6),
            facecolors=["#FFFFFF"],
            edgecolors=["#000000"],
        )
        if each_chrom_df.empty:
            continue

        plot_x = list(
            each_chrom_df.apply(lambda x: (int(x.start), int(x.span)), axis=1)
        )
        # plot_colors = list(each_chrom_df["Origin"].map(lambda x: COLOR_MAP[x]))
        plot_colors = ["#43CD80"] * len(each_chrom_df)
        ax.broken_barh(plot_x, (y_pos, 6), facecolors=plot_colors)
    ax.set_yticks(
        [i * 10 + 7 for i in range(plot_height - 1)],
        labels=np.flip(chr_df.index),
    )
    ax.set_xticks(
        xaxis_ticks,
        xaxis_labels,
    )

    out_png = out_dir / f"{name}.png"
    out_pdf = out_png.with_suffix(".pdf")
    plot_title = f"{name}"
    plt.title(plot_title)
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()


def main(
    input_dir: Path,
    chrom_sizes: Path,
    output_dir: Path,
    min_coverage: int = 10,
    split_bed: Optional[Path] = None,
    sample_list: Optional[Path] = None,
) -> None:
    """
    Plot coverage plot for each sample in input_dir.

    Parameters
    ----------
    input_dir : Path
        Path to directory contains coverage bed files.
    chrom_sizes : Path
        Path to chromosome size file.
    output_dir : Path
        Path to output directory.
    min_coverage : int, default=10
        Minimum coverage to include in plot.
    split_bed : Path, optional
        Path to split_bed file. If not None, merge chromosomes using split_bed
        file.
    sample_list : Path, optional
        Path to sample list file. If not None, only process samples in sample
        list file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    bed_files = list(input_dir.glob("*.bed"))
    chrom_df = pd.read_table(
        chrom_sizes, header=None, names=["chrom", "length"], index_col=0
    )
    for bed_file in tqdm(bed_files):
        sample_name = bed_file.stem.replace(".cov", "")
        if sample_list and sample_name not in list(
            pd.read_table(sample_list, header=None)[0]
        ):
            continue
        cov_df = pd.read_table(
            bed_file, header=None, names=["chrom", "start", "end", "cov"]
        )
        cov_df["length"] = cov_df["end"] - cov_df["start"]
        cov_df = cov_df[cov_df["cov"] >= min_coverage * cov_df["length"]]
        if split_bed:
            cov_df = merge_chr(cov_df, split_bed)
        plot_cov(cov_df, chrom_df, output_dir, sample_name)


if __name__ == "__main__":
    typer.run(main)
