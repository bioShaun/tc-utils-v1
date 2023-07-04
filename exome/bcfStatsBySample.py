import typer
import pandas as pd

from io import StringIO
from pathlib import Path

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
    "nRefHom",
    "nNonRefHom",
    "nHets",
    "nIndels",
    "nMissing",
    "nTransitions",
    "nTransversions",
    "averageDepth",
]


def sampleStats(bcfstats: Path, out_file: Path) -> None:
    bcfstats_list = open(bcfstats).readlines()
    psc_stats = "".join([each for each in bcfstats_list if each.startswith("PSC")])
    psc_df = pd.read_csv(StringIO(psc_stats), sep="\t", header=None, names=ALL_COLUMNS)
    psc_df.to_csv(out_file, sep="\t", index=False, columns=OUT_COLUMNS)


if __name__ == "__main__":
    typer.run(sampleStats)
