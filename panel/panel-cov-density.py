from enum import Enum
from pathlib import Path
from typing import List, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer
from loguru import logger
from tqdm import tqdm


class ColorMap(str, Enum):
    RdYlGn = "RdYlGn_r"
    RdYlBu = "RdYlBu_r"
    Greens = "Greens"
    Reds = "Reds"
    YlGn = "YlGn"


def add_region(
    df: pd.DataFrame,
    chr_df: pd.DataFrame,
    region_list: List[int],
) -> pd.DataFrame:
    add_region_df_list = []
    for chrom, each_chrom_df in df.groupby("chrom"):
        logger.info(f"add region for chrom {chrom}")
        # if chrom not in chr_df.index:
        # continue
        chrom_len = chr_df.loc[chrom, "chrom_length"]  # type: ignore
        for region_size in region_list:
            region_size_bp = region_size * 1000
            each_chrom_cut = pd.cut(
                each_chrom_df["pos"],
                range(0, chrom_len + region_size_bp, region_size_bp),
            )
            each_chrom_df[f"region_{region_size}k"] = each_chrom_cut.map(
                lambda x: f"{chrom}_{x.left}_{x.right}"
            )
        add_region_df_list.append(each_chrom_df)
    if add_region_df_list:
        return pd.concat(add_region_df_list)
    return pd.DataFrame([])


def make_window(chr_size_df: pd.DataFrame, region_size: int) -> pd.DataFrame:
    df_list = []
    region_size = region_size * 1000
    for row in chr_size_df.itertuples():
        for i in range(0, row.chrom_length, region_size):
            df_list.append([row.Index, i, i + region_size])
    return pd.DataFrame(df_list, columns=["chrom", "start", "end"])


def cal_region_coverage(probe_df: pd.DataFrame, chr_df: pd.DataFrame, region_size: int):
    plot_df = probe_df[probe_df["chrom"].isin(chr_df.index)]
    region_label = f"region_{region_size}k"
    region_count_df = pd.DataFrame(plot_df.groupby(["chrom", region_label]).size())
    region_count_df.columns = ["coverage"]
    region_count_df.reset_index(inplace=True)
    region_count_df["end"] = region_count_df[region_label].map(
        lambda x: x.split("_")[-1]
    )
    region_count_df["start"] = region_count_df[region_label].map(
        lambda x: x.split("_")[-2]
    )
    region_count_df[["start", "end"]] = region_count_df[["start", "end"]].astype(int)
    region_count_df = region_count_df[["chrom", "start", "end", "coverage"]]
    region_count_df["coverage"] = region_count_df["coverage"].astype("int")
    window_df = make_window(chr_df, region_size)
    region_count_df = window_df.merge(region_count_df, how="left").fillna(0)
    region_count_df["coverage"] = region_count_df["coverage"].astype("int")
    return region_count_df


def probe_count_tag(count: int) -> str:
    if count < 1000:
        return f"{count}K"
    else:
        return f"{count / 1000:.1f}M"


def get_coverage_rate(region_count_df: pd.DataFrame) -> float:
    covered_df = region_count_df[region_count_df["coverage"] > 0]
    coverage_rate = round(100 * len(covered_df) / len(region_count_df), 1)
    return coverage_rate


def get_plot_size(df: pd.DataFrame) -> Tuple[float, int]:
    width = max(df.groupby("chrom").size().max() * 0.0075, 24)  # type: ignore
    height = df["chrom"].unique().size * 1
    if width > height * 2:
        width = height * 2
    return width, height


def get_plot_xaxis(df: pd.DataFrame) -> Tuple[List[int], List[str]]:
    chrom_max_length = df["end"].max()
    mega_base = np.floor(np.log10(chrom_max_length))
    max_chr_show_length = np.ceil(chrom_max_length / 10**mega_base)
    x_axis_ticks = [int(i * 10**mega_base) for i in range(int(max_chr_show_length) + 1)]
    x_axis_labels = [
        f"{int(i * 10**mega_base / 1e6)}M" for i in range(int(max_chr_show_length) + 1)
    ]
    return x_axis_ticks, x_axis_labels


def density_plot(
    density_df: pd.DataFrame,
    density_path: Path,
    plot_title: str,
    max_coverage: int,
    color_map: ColorMap = ColorMap.RdYlGn,
    min_coverage: int = 1,
) -> None:
    plt.rcParams["font.size"] = 24
    cmap = mpl.colormaps[color_map].with_extremes(over="0.25", under="0.75")  # type: ignore
    norm = mpl.colors.Normalize(vmin=min_coverage, vmax=max_coverage)
    fig, ax = plt.subplots()
    plot_width, plot_height = get_plot_size(density_df)
    fig.set_figheight(plot_height)
    fig.set_figwidth(plot_width)

    chroms = density_df["chrom"].unique()

    for n, chrom in enumerate(chroms):
        df = density_df[density_df["chrom"] == chrom]
        plot_data = norm(df["coverage"].values).data
        plot_x = list(df.apply(lambda x: (x.start, x.end - x.start), axis=1))
        y_pos = (plot_height - n - 1) * 10 + 4
        ax.broken_barh(plot_x, (y_pos, 6), facecolors=cmap(plot_data))
    ax.set_yticks(
        [i * 10 + 7 for i in range(plot_height)],
        labels=np.flip(chroms),
    )
    xaxis_ticks, xaxis_labels = get_plot_xaxis(density_df)
    ax.set_xticks(
        xaxis_ticks,
        xaxis_labels,
    )
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)  # type: ignore
    sm.set_array([])
    plt.title(plot_title)
    fig.colorbar(sm, ax=ax, label="")
    plt.savefig(density_path.with_suffix(".density.png"), dpi=300)
    plt.savefig(density_path.with_suffix(".density.pdf"))


def plot_probe_coverage(
    probe_df: pd.DataFrame,
    region_size: int,
    chr_size_df: pd.DataFrame,
    out_dir: Path,
    probe_name: str,
    color_map: ColorMap,
    max_coverage: int,
    min_coverage: int = 0,
) -> None:
    region_count_df = cal_region_coverage(
        probe_df=probe_df, chr_df=chr_size_df, region_size=region_size
    )
    region_tag = probe_count_tag(region_size)
    region_count_file = out_dir / f"{probe_name}.{region_tag}-window.coverage.tsv"
    region_count_df.to_csv(region_count_file, sep="\t", index=False)
    plot_title = f"{probe_name} ({region_tag} window coverage)"
    density_plot(
        density_df=region_count_df,
        density_path=region_count_file,
        plot_title=plot_title,
        color_map=color_map,
        max_coverage=max_coverage,
        min_coverage=min_coverage,
    )


def main(
    cov_file_dir: Path,
    chr_size_file: Path,
    region_size: int,
    region_count: int,
    out_dir: Path,
    color_map: ColorMap = ColorMap.RdYlGn,
    coverage_cutoff: int = 10,
    min_coverage: int = 0,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    bed_files = list(cov_file_dir.glob("*.bed"))
    chr_size_df = pd.read_table(
        chr_size_file,
        header=None,
        names=["chrom", "chrom_length"],
        index_col=0,
        usecols=[0, 1],
    )
    for bed_file in tqdm(bed_files):
        sample_name = bed_file.stem.replace(".cov", "")
        probe_df = pd.read_table(
            bed_file,
            header=None,
            usecols=[0, 1, 2, 3],
            names=["chrom", "start", "end", "coverage"],
        )
        probe_df["pos"] = (probe_df["start"] + probe_df["end"]) // 2
        probe_df["span"] = probe_df["end"] - probe_df["start"]

        probe_df["chrom"] = probe_df["chrom"].astype(str)
        chr_size_df.index = chr_size_df.index.astype(str)
        probe_df = probe_df[probe_df["chrom"].isin(chr_size_df.index)]
        probe_df = probe_df[
            probe_df["coverage"] >= coverage_cutoff * probe_df["span"]
        ].copy()
        probe_df = add_region(probe_df, chr_size_df, [region_size])
        plot_probe_coverage(
            probe_df=probe_df,
            chr_size_df=chr_size_df,
            region_size=region_size,
            out_dir=out_dir,
            probe_name=sample_name,
            color_map=color_map,
            min_coverage=min_coverage,
            max_coverage=region_count,
        )


if __name__ == "__main__":
    typer.run(main)
