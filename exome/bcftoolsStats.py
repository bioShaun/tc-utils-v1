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


INDEL_COLUMNS = [
    "PSI",
    "id",
    "Sample",
    "in_frame",
    "out_frame",
    "not_applicable",
    "out_vs_in_out_ratio",
    "nInsHets",
    "nDelHets",
    "nInsAltHoms",
    "nDelAltHoms",
]

OUT_COLUMNS = [
    "Sample",
    "nRefHom",
    "nNonRefHom",
    "nHets",
    "nIndels",
    "nMissing",
    "nTransitions",
    "nTransversions",
    "averageDepth",
]

INDEL_OUT_COLUMNS = [
    "Sample",
    "nRefHom",
    "nInsHets",
    "nDelHets",
    "nInsAltHoms",
    "nDelAltHoms",
    "nMissing",
    "averageDepth",
]


def sampleStats(bcfstats: Path, out_file: Path, indel: bool = False) -> None:
    bcfstats_list = open(bcfstats).readlines()
    psc_stats = "".join([each for each in bcfstats_list if each.startswith("PSC")])
    psc_df = pd.read_csv(StringIO(psc_stats), sep="\t", header=None, names=ALL_COLUMNS)
    out_cols = OUT_COLUMNS[:]
    if indel:
        psi_stats = "".join([each for each in bcfstats_list if each.startswith("PSI")])
        psi_df = pd.read_csv(
            StringIO(psi_stats), sep="\t", header=None, names=INDEL_COLUMNS
        )
        psi_df.drop(["PSI", "id"], axis=1, inplace=True)
        out_cols = INDEL_COLUMNS[:]
        psc_df = psc_df.merge(psi_df, on="Sample", how="left")
        psc_df["Total"] = psc_df[
            [
                "nRefHom",
                "nInsHets",
                "nDelHets",
                "nInsAltHoms",
                "nDelAltHoms",
                "nMissing",
            ]
        ].sum(1)
    psc_df = psc_df[out_cols]
    psc_df.to_csv(out_file, sep="\t", index=False, columns=OUT_COLUMNS)


if __name__ == "__main__":
    typer.run(sampleStats)
