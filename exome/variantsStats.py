from io import StringIO
from pathlib import Path

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

SNP_OUT_COLUMNS = [
    "Total",
    "nRefHom",
    "nNonRefHom",
    "nHets",
    "nMissing",
    "nSingletons",
    "nTransitions",
    "nTransversions",
]

PSI_COLUMNS = [
    "PSI",
    "id",
    "Sample",
    "inFrame",
    "outFrame",
    "notApplicable",
    "outRatio",
    "nInsHets",
    "nDelHets",
    "nInsAltHoms",
    "nDelAltHoms",
]

INDEL_OUT_COLUMNS = [
    "Total",
    "nRefHom",
    "nInsHets",
    "nDelHets",
    "nInsAltHoms",
    "nDelAltHoms",
    "nMissing",
    "nSingletons",
]


def snpStats(bcfstats: Path) -> pd.DataFrame:
    bcfstats_list = open(bcfstats).readlines()
    psc_stats = "".join([each for each in bcfstats_list if each.startswith("PSC")])
    psc_df = pd.read_csv(StringIO(psc_stats), sep="\t", header=None, names=ALL_COLUMNS)
    psc_df["Total"] = psc_df[["nRefHom", "nNonRefHom", "nHets", "nMissing"]].sum(1)
    psc_df = psc_df.set_index("Sample")
    psc_df = psc_df[SNP_OUT_COLUMNS].copy()
    multi_columns = [["SNP"] * len(psc_df.columns), list(psc_df.columns)]
    psc_df.columns = pd.MultiIndex.from_arrays(multi_columns)
    return psc_df


def indelStats(indel_stats: Path) -> pd.DataFrame:
    bcfstats_list = open(indel_stats).readlines()
    psc_stats = "".join([each for each in bcfstats_list if each.startswith("PSC")])
    psc_df = pd.read_csv(StringIO(psc_stats), sep="\t", header=None, names=ALL_COLUMNS)
    psc_df["Total"] = psc_df[["nRefHom", "nIndels", "nMissing"]].sum(1)
    psc_df.drop(columns=["id"], inplace=True)
    psi_stats = "".join([each for each in bcfstats_list if each.startswith("PSI")])
    psi_df = pd.read_csv(StringIO(psi_stats), sep="\t", header=None, names=PSI_COLUMNS)
    psi_df.drop(columns=["id"], inplace=True)
    merged_df = pd.merge(psc_df, psi_df, on=["Sample"])
    merged_df = merged_df.set_index("Sample")
    merged_df = merged_df[INDEL_OUT_COLUMNS].copy()
    multi_columns = [["INDEL"] * len(merged_df.columns), list(merged_df.columns)]
    merged_df.columns = pd.MultiIndex.from_arrays(multi_columns)
    return merged_df


def sampleStats(snp_stats: Path, indel_stats: Path, out_file: Path):
    snp_df = snpStats(snp_stats)
    indel_df = indelStats(indel_stats)
    merged_df = pd.merge(snp_df, indel_df, left_index=True, right_index=True)
    merged_df.reset_index(col_level=1, inplace=True)
    merged_df.to_csv(
        out_file,
    )


if __name__ == "__main__":
    typer.run(sampleStats)
