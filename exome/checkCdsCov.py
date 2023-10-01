import typer
from pathlib import Path
import pandas as pd
from functools import reduce
from typing import List, Optional, Tuple
from loguru import logger
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from tqdm import tqdm

BED_COLUMNS = ["chrom", "start", "end", "transcript"]
# BED_COLUMNS = ["chrom", "start", "end"]


def merge_chr(df: pd.DataFrame, split_bed: Path) -> pd.DataFrame:
    split_bed_df = pd.read_csv(
        split_bed,
        header=None,
        names=["new_chrom", "offset", "offset_end", "chrom"],
        sep="\t",
    )
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


def load_bed_files(
    bed_dir: Path, bed_cols: List[str], split_bed: Optional[Path] = None
) -> pd.DataFrame:
    df_list = []
    bed_files = list(bed_dir.glob("*.bed"))
    for bed_i in tqdm(bed_files):
        # logger.info(f"Load {bed_i} ...")
        sample_name = bed_i.stem.rstrip(".cov")
        bed_i_columns = [*bed_cols, sample_name]
        df_i = pd.read_table(bed_i, header=None, names=bed_i_columns)
        df_list.append(df_i)
    df = reduce(
        lambda x, y: pd.merge(x, y),
        df_list,
    )
    if split_bed is None:
        return df
    else:
        return merge_chr(df, split_bed)


def add_region(df: pd.DataFrame, chr_size: Path, region_size: int) -> pd.DataFrame:
    chr_df = pd.read_table(
        chr_size,
        header=None,
        names=["chrom", "chrom_length"],
        index_col=0,
    )
    add_region_df_list = []
    for chrom, each_chrom_df in df.groupby("chrom"):
        logger.info(f"add region for chrom {chrom}")
        if chrom not in chr_df.index:
            continue
        chrom_len = chr_df.loc[chrom, "chrom_length"]  # type: ignore
        each_chrom_cut = pd.cut(
            each_chrom_df["pos"],
            range(0, chrom_len + region_size, region_size),
        )
        each_chrom_df["region"] = each_chrom_cut.map(
            lambda x: f"{chrom}_{x.left}_{x.right}"
        )
        add_region_df_list.append(each_chrom_df)
    if add_region_df_list:
        return pd.concat(add_region_df_list)
    return pd.DataFrame([])


def get_coverage_rate(region_count_df: pd.DataFrame) -> float:
    miss_df = region_count_df[region_count_df["miss"] == 1]
    coverage_rate = round(100 * len(miss_df) / len(region_count_df), 1)
    return coverage_rate


def get_plot_size(df: pd.DataFrame) -> Tuple[float, int]:
    width = df.groupby("chrom").size().max() * 0.0075
    height = df["chrom"].unique().size * 1
    if width < 24:
        width = 24
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


def density_plot(
    density_df: pd.DataFrame,
    density_path: Path,
    plot_title: str,
    plot_type: str,
) -> None:
    plt.rcParams["font.size"] = 24
    # .with_extremes(over="0.25", under="0.75")  # type: ignore
    # norm = mpl.colors.Normalize(vmin=1, vmax=density_df["coverage"].max())
    if plot_type == "miss":
        plot_lable = "miss"
        color_map = "Reds"
        cmap = ListedColormap(["silver", "red"])
    else:
        plot_lable = "coverage"
        color_map = "RdYlGn"
        cmap = mpl.colormaps[color_map]
    fig, ax = plt.subplots()
    plot_width, plot_height = get_plot_size(density_df)
    fig.set_figheight(plot_height)
    fig.set_figwidth(plot_width)

    chroms = density_df["chrom"].unique()

    for n, chrom in enumerate(chroms):
        df = density_df[density_df["chrom"] == chrom]
        plot_data = df[plot_lable].values
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
    plt.title(plot_title)
    if plot_type != "miss":
        sm = plt.cm.ScalarMappable(cmap=cmap)  # type: ignore
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label="")
    plt.savefig(f"{density_path}.density.png", dpi=300)
    plt.savefig(f"{density_path}.density.pdf")


def plot_probe_coverage(
    df: pd.DataFrame,
    region_size: int,
    out_dir: Path,
    plot_type="miss",
    sample_name: Optional[str] = None,
) -> None:
    df[["chrom", "start", "end"]] = df["region"].str.split("_", expand=True)
    df[["start", "end"]] = df[["start", "end"]].astype("int")
    df.drop("region", inplace=True, axis=1)

    region_tag = region_size // 1_000_000
    if plot_type == "miss":
        miss_rate = get_coverage_rate(df)
        plot_title = f"{region_tag}M window miss rate: {miss_rate}%"
    else:
        plot_title = f"{region_tag}M Coverage"
    if plot_type == "miss":
        out_name = f"{region_tag}M-miss"
    else:
        out_name = f"{region_tag}M-coverage"
    if sample_name is not None:
        out_name = f"{sample_name}-{out_name}"
        plot_title = f"{sample_name} {plot_title}"
    out_prefix = out_dir / out_name
    df.to_csv(
        f"{out_prefix}.tsv",
        sep="\t",
        index=False,
        columns=["chrom", "start", "end", plot_type],
    )
    density_plot(
        density_df=df,
        density_path=out_prefix,
        plot_title=plot_title,
        plot_type=plot_type,
    )


def main(
    cds_cov_dir: Path,
    out_file: Path,
    min_reads: int = 10,
    split_bed: Optional[Path] = None,
    transcript_id: bool = False,
) -> None:
    if transcript_id:
        loci_columns = BED_COLUMNS[:]
    else:
        loci_columns = BED_COLUMNS[:-1]
    cds_df = load_bed_files(cds_cov_dir, loci_columns, split_bed)
    df_matrix = cds_df.set_index(loci_columns)
    df_matrix_bool = df_matrix >= min_reads
    cov_ratio = df_matrix_bool.sum() / df_matrix.shape[0]
    cov_ratio.to_csv(out_file, sep="\t", header=False, float_format="%.3f")


if __name__ == "__main__":
    typer.run(main)
