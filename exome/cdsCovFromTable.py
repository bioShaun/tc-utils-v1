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

BED_COLUMNS = ["chrom", "start", "end", "transcript"]


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


def load_bed_files(bed_dir: Path, split_bed: Optional[Path] = None) -> pd.DataFrame:
    df_list = []
    for bed_i in bed_dir.glob("*.bed"):
        logger.info(f"Load {bed_i} ...")
        sample_name = bed_i.stem.rstrip(".cov")
        bed_i_columns = [*BED_COLUMNS, sample_name]
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
    cds_cov_table: Path,
    sample_file: Path,
    out_dir: Path,
    chr_size: Path,
    min_reads: int = 10,
    region_size: int = 1000000,
    sample_ratio_cutoff: float = 0.5,
    region_ratio_cutoff: float = 0.5,
    split_bed: Optional[Path] = None,
    plot_all: bool = False,
    transcript_id: bool = False,
) -> None:
    out_dir.mkdir(exist_ok=True, parents=True)
    if transcript_id:
        loci_columns = BED_COLUMNS[:]
    else:
        loci_columns = BED_COLUMNS[:-1]
    sample_list = [each.strip() for each in sample_file.open()]
    cds_columns = [*loci_columns, *sample_list]
    cds_df = pd.read_table(cds_cov_table, header=None, names=cds_columns)
    if split_bed is not None:
        cds_df = merge_chr(cds_df, split_bed)
    df_matrix = cds_df.set_index(loci_columns)
    df_matrix_bool = df_matrix >= min_reads
    cover_df = df_matrix_bool.sum(1)
    cover_ratio_df = cover_df / df_matrix_bool.shape[1]
    passed_df_matrix_bool = df_matrix_bool.loc[
        cover_ratio_df[cover_ratio_df >= sample_ratio_cutoff].index
    ].copy()
    if plot_all:
        cover_ratio_df.name = "cover_ratio"
        cover_ratio_df = pd.DataFrame(cover_ratio_df).reset_index()
        cover_ratio_df["pos"] = (cover_ratio_df["start"] + cover_ratio_df["end"]) // 2
        cover_ratio_df["cover_passed"] = (
            cover_ratio_df["cover_ratio"] >= sample_ratio_cutoff
        )
        cover_ratio_df = add_region(cover_ratio_df, chr_size, region_size)
        region_passed = cover_ratio_df.groupby("region")["cover_passed"].sum()
        region_cds_count = cover_ratio_df.groupby("region").size()
        region_ratio = region_passed / region_cds_count
        region_ratio.dropna(inplace=True)
        region_ratio.name = "coverage"
        region_miss = region_ratio.map(lambda x: 0 if x >= region_ratio_cutoff else 1)
        region_miss.name = "miss"
        region_ratio_df = pd.DataFrame(region_ratio).reset_index()
        region_miss_df = pd.DataFrame(region_miss).reset_index()
        plot_probe_coverage(df=region_miss_df, region_size=region_size, out_dir=out_dir)
        plot_probe_coverage(
            df=region_ratio_df,
            region_size=region_size,
            out_dir=out_dir,
            plot_type="coverage",
        )
    for sample_i in passed_df_matrix_bool.columns:
        logger.info(f"Processing {sample_i}")
        sampe_i_df = passed_df_matrix_bool[[sample_i]].reset_index()
        sampe_i_df["pos"] = (sampe_i_df["start"] + sampe_i_df["end"]) // 2
        cover_ratio_df = add_region(sampe_i_df, chr_size, region_size)
        region_passed = cover_ratio_df.groupby("region")[sample_i].sum()
        region_cds_count = cover_ratio_df.groupby("region").size()
        region_ratio = region_passed / region_cds_count
        region_ratio.dropna(inplace=True)
        region_ratio.name = "coverage"
        region_miss = region_ratio.map(lambda x: 0 if x >= region_ratio_cutoff else 1)
        region_miss.name = "miss"
        region_ratio_df = pd.DataFrame(region_ratio).reset_index()
        region_miss_df = pd.DataFrame(region_miss).reset_index()
        plot_probe_coverage(
            df=region_miss_df,
            region_size=region_size,
            out_dir=out_dir,
            sample_name=sample_i,
        )
        plot_probe_coverage(
            df=region_ratio_df,
            region_size=region_size,
            out_dir=out_dir,
            plot_type="coverage",
            sample_name=sample_i,
        )


if __name__ == "__main__":
    typer.run(main)
