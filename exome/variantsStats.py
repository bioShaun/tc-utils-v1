from io import StringIO
from pathlib import Path
from typing import List

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


def extract_coverage_stats(alignment_stats: Path, region_size: int) -> List[float]:
    out_stats = [0.0, 0.0]
    with alignment_stats.open() as align_info:
        for each_line in align_info:
            if "bases mapped (cigar)" in each_line:
                alignment_length = int(each_line.split("\t")[2])
                out_stats[0] = alignment_length / region_size
            if "percentage of target genome with coverage" in each_line:
                coverage = float(each_line.split("\t")[2])
                out_stats[1] = coverage
    return out_stats


def coverageStats(alignment_dir: Path, region_size: int) -> pd.DataFrame:
    cov_data = []
    for cov_file in alignment_dir.glob("*/*.coverage.stat"):
        sample_name = cov_file.parent.name
        cov_percent, cov_depth = extract_coverage_stats(cov_file, region_size)
        cov_data.append([sample_name, cov_percent, cov_depth])
    cov_df = pd.DataFrame(
        cov_data, columns=["Sample", "Mean Sequencing Depth", "Coverage >= 1x(%)"]
    )
    cov_df = cov_df.set_index("Sample")
    multi_columns = [[""] * len(cov_df.columns), list(cov_df.columns)]
    cov_df.columns = pd.MultiIndex.from_arrays(multi_columns)
    return cov_df


def get_region_size(bed_file: Path) -> int:
    bed_df = pd.read_table(bed_file, header=None, names=["chrom", "start", "end"])
    bed_df["region_size"] = bed_df["end"] - bed_df["start"]
    return bed_df["region_size"].sum()


def sampleStats(snp_stats: Path, indel_stats: Path, region_bed: Path, out_file: Path):
    snp_df = snpStats(snp_stats)
    indel_df = indelStats(indel_stats)
    region_size = get_region_size(region_bed)
    cov_df = coverageStats(snp_stats.parent, region_size)
    merged_df = pd.merge(snp_df, indel_df, left_index=True, right_index=True)
    merged_df = pd.merge(cov_df, merged_df, left_index=True, right_index=True)
    merged_df.reset_index(col_level=1, inplace=True)
    merged_df.to_csv(out_file, index=False)


if __name__ == "__main__":
    typer.run(sampleStats)
