from pathlib import Path
from typing import List, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer


def get_plot_size(df: pd.DataFrame) -> Tuple[float, int]:
    width = df.groupby("chrom").size().max() * 0.0125
    if width < 15:
        width = 15
    height = df["chrom"].unique().size * 1
    return width, height


def get_plot_xaxis(df: pd.DataFrame) -> Tuple[List[int], List[str]]:
    chrom_max_length = df["end"].max()
    mega_base = np.floor(np.log10(chrom_max_length))
    max_chr_show_length = np.ceil(chrom_max_length / 10**mega_base)
    x_axis_ticks = [
        int(i * 10**mega_base) for i in range(int(max_chr_show_length) + 1)
    ]
    x_axis_labels = [
        f"{int(i * 10**mega_base / 1e6)}M" for i in range(int(max_chr_show_length) + 1)
    ]
    return x_axis_ticks, x_axis_labels


def main(density_file: Path, plot_title: str, log_scale: bool = False) -> None:
    plt.rcParams["font.size"] = 24

    density_df = pd.read_table(
        density_file, header=None, names=["chrom", "start", "end", "density"]
    )
    cmap = mpl.colormaps["RdYlGn_r"].with_extremes(over="0.25", under="0.75")
    if log_scale:
        density_df["density"] = density_df["density"].map(lambda x: np.log2(x + 1))
    norm = mpl.colors.Normalize(vmin=1, vmax=density_df["density"].max())
    fig, ax = plt.subplots()
    plot_width, plot_height = get_plot_size(density_df)
    fig.set_figheight(plot_height)
    fig.set_figwidth(plot_width)

    for n, (_, df) in enumerate(density_df.groupby("chrom")):
        plot_data = norm(df["density"].values).data
        plot_x = list(df.apply(lambda x: (x.start, x.end - x.start), axis=1))
        y_pos = (plot_height - n - 1) * 10 + 4
        ax.broken_barh(plot_x, (y_pos, 6), facecolors=cmap(plot_data))
    ax.set_yticks(
        [i * 10 + 7 for i in range(plot_height)],
        labels=np.flip(density_df["chrom"].unique()),
    )
    xaxis_ticks, xaxis_labels = get_plot_xaxis(density_df)
    ax.set_xticks(
        xaxis_ticks,
        xaxis_labels,
    )
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.title(plot_title)
    fig.colorbar(sm, ax=ax, label="")
    plt.savefig(density_file.with_suffix(".png"), dpi=300)
    plt.savefig(density_file.with_suffix(".pdf"))


if __name__ == "__main__":
    typer.run(main)
