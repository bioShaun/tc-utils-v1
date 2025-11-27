from io import StringIO
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import typer

ALL_COLUMNS = [
    "PSC",
    "id",
    "Sample",
    "nRefHom",
    "nNonRefHom",
    "nHets",
    "nTransitions",
    "nTransversions",
    "nIndels",
    "averageDepth",
    "nSingletons",
    "nHapRef",
    "nHapAlt",
    "nMissing",
]

OUT_COLUMNS = [
    "Sample",
    "NON-MISS",
    "REF",
    "ALT",
    "HET",
    "MISS_RATIO",
    "REF_RATIO",
    "ALT_RATIO",
    "HET_RATIO",
]


COLUMN_MAP = {
    "nRefHom": "REF",
    "nNonRefHom": "ALT",
    "nHets": "HET",
}


def plot(df: pd.DataFrame, outdir: Path) -> None:
    # 第一张图：MISS_RATIO
    df.sort_values(by="MISS_RATIO", inplace=True, ascending=False)
    df = df.reset_index(drop=True)
    df["Index"] = df.index + 1
    plt.figure(figsize=(12, 5))
    plt.plot(df["Index"], df["MISS_RATIO"], color="red")
    # plt.xlabel("Sample Index")
    plt.ylabel("MISS_RATIO")
    plt.title("MISS_RATIO across Samples")
    # plt.xticks(range(0, len(df) + 1, 50))

    # y轴：最高 0.1 或最大值
    ymax = max(0.25, df["MISS_RATIO"].max())
    plt.ylim(0, ymax)

    miss_plot_file_pdf = outdir / "MISS_RATIO.pdf"
    miss_plot_file_png = outdir / "MISS_RATIO.png"
    plt.tight_layout()
    plt.savefig(miss_plot_file_pdf)
    plt.savefig(miss_plot_file_png, dpi=300)
    plt.close()

    # 第二张图：HET_RATIO
    df.sort_values(by="HET_RATIO", inplace=True, ascending=False)
    df = df.reset_index(drop=True)
    df["Index"] = df.index + 1
    het_plot_file_pdf = outdir / "HET_RATIO.pdf"
    het_plot_file_png = outdir / "HET_RATIO.png"
    plt.figure(figsize=(12, 5))
    plt.plot(df["Index"], df["HET_RATIO"], color="blue")
    # plt.xlabel("Sample Index")
    plt.ylabel("HET_RATIO")
    plt.title("HET_RATIO across Samples")
    # plt.xticks(range(0, len(df) + 1, 50))

    # y轴固定到 1
    plt.ylim(0, 1)

    plt.tight_layout()
    plt.savefig(het_plot_file_pdf)
    plt.savefig(het_plot_file_png, dpi=300)
    plt.close()


def sampleStats(bcfstats: Path, out_file: Path) -> None:
    out_dir = out_file.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    bcfstats_list = open(bcfstats).readlines()
    psc_stats = "".join([each for each in bcfstats_list if each.startswith("PSC")])
    psc_df = pd.read_csv(StringIO(psc_stats), sep="\t", header=None, names=ALL_COLUMNS)
    psc_df["Total"] = (
        psc_df["nRefHom"] + psc_df["nNonRefHom"] + psc_df["nHets"] + psc_df["nMissing"]
    )
    psc_df["NON-MISS"] = psc_df["Total"] - psc_df["nMissing"]
    psc_df.rename(columns=COLUMN_MAP, inplace=True)
    psc_df["MISS_RATIO"] = psc_df["nMissing"] / psc_df["Total"]
    psc_df["REF_RATIO"] = psc_df["REF"] / psc_df["Total"]
    psc_df["ALT_RATIO"] = psc_df["ALT"] / psc_df["Total"]
    psc_df["HET_RATIO"] = psc_df["HET"] / psc_df["Total"]
    psc_df.to_excel(out_file, index=False, columns=OUT_COLUMNS, float_format="%.4f")
    plot(psc_df, out_file.parent)


if __name__ == "__main__":
    typer.run(sampleStats)
