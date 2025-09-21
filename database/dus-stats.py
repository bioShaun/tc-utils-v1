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


def sampleStats(bcfstats: Path, out_file: Path) -> None:
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
    psc_df.to_csv(out_file, sep="\t", index=False, columns=OUT_COLUMNS)


if __name__ == "__main__":
    typer.run(sampleStats)
