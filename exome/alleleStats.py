import typer

import numpy as np
import pandas as pd

from pathlib import Path


app = typer.Typer()


def hetRatio(obs_stats: str) -> float:
    hom1, het, hom2 = [float(each) for each in obs_stats.split("/")]
    return het / (hom1 + het + hom2)


def ref_freq(row) -> float:
    hom1 = row.value_counts().get("0/0", 0)
    het = row.value_counts().get("0/1", 0)
    hom2 = row.value_counts().get("1/1", 0)
    return (2 * hom1 + het) / (2 * (hom1 + het + hom2))


def het_freq(row) -> float:
    hom1 = row.value_counts().get("0/0", 0)
    het = row.value_counts().get("0/1", 0)
    hom2 = row.value_counts().get("1/1", 0)
    return het / (hom1 + het + hom2)


def shannonIndex(p_ref: float, p_alt: float) -> float:
    if p_ref == 0 or p_alt == 0:
        return 0
    else:
        return -1 * (p_ref * np.log(p_ref) + p_alt * np.log(p_alt))


@app.command()
def by_locus(vcftools_frq: Path, vcftools_hwe: Path, out_file: Path) -> None:
    df = pd.read_csv(vcftools_frq, sep="\t")
    df = df.reset_index()
    df.columns = ["CHR", "POS", "N_ALLELES", "N_CHR", "REF", "ALT"]
    hwe_df = pd.read_csv(vcftools_hwe, sep="\t", usecols=[0, 1, 2, 3])
    df = hwe_df.merge(df)
    df["maf"] = df.apply(lambda row: max(row["REF"], row["ALT"]), axis=1)
    df["Ho"] = df.apply(lambda row: hetRatio(row["OBS(HOM1/HET/HOM2)"]), axis=1)
    df["He"] = df.apply(lambda row: hetRatio(row["E(HOM1/HET/HOM2)"]), axis=1)
    df["PIC"] = df["maf"].map(
        lambda maf: 1 - (maf**2 + (1 - maf) ** 2) - 2 * (maf**2) * ((1 - maf) ** 2)
    )
    df["SHI"] = df.apply(lambda row: shannonIndex(row["REF"], row["ALT"]), axis=1)
    df.to_csv(out_file, index=False)


@app.command()
def by_sample(gt_tsv: Path, sample_list: Path, out_file: Path) -> None:
    df = pd.read_csv(gt_tsv, sep="\t", header=None)
    df = df.drop(columns=[0, 1, 2, 3])
    df.columns = pd.read_csv(sample_list, sep="\t", header=None)[0].tolist()
    df = df.T
    df["p1"] = df.apply(lambda row: ref_freq(row), axis=1)
    df["p2"] = 1 - df["p1"]
    df["Ho"] = df.apply(lambda row: het_freq(row), axis=1)
    df["He"] = df.apply(lambda row: 1 - row["p1"] ** 2 - row["p2"] ** 2, axis=1)
    df["PIC"] = df.apply(
        lambda row: 1
        - (row["p1"] ** 2 + row["p2"] ** 2)
        - 2 * (row["p1"] ** 2) * (row["p2"] ** 2),
        axis=1,
    )
    df["SHI"] = df.apply(lambda row: shannonIndex(row["p1"], row["p2"]), axis=1)
    df.to_csv(out_file, columns=["p1", "p2", "Ho", "He", "PIC", "SHI"])


if __name__ == "__main__":
    app()
