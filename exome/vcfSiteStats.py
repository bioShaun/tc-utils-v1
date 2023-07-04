from enum import Enum
from functools import partial
from pathlib import Path

import pandas as pd
import typer


class StatType(str, Enum):
    HET = "het"
    MISS = "miss"
    MAF = "maf"


def allele_stats(row: pd.Series, stat_type: StatType) -> float:
    allele_count = row.value_counts()
    miss_count = allele_count.get("./.", 0)
    het_count = allele_count.get("0/1", 0)
    real_sample_count = len(row) - miss_count
    match stat_type:
        case StatType.MISS:
            return miss_count / len(row)
        case StatType.HET:
            return het_count / real_sample_count
        case StatType.MAF:
            ref_count = allele_count.get("0/0", 0)
            ref_rate = (ref_count * 2 + het_count) / (real_sample_count * 2)
            maf = ref_rate if ref_rate < 0.5 else 1 - ref_rate
            return maf
        case _:
            raise ValueError(f"Unknown type {stat_type}")


def transform_one(df: pd.DataFrame) -> pd.DataFrame:
    melt_df = df.melt(id_vars=["CHROM", "POS", "REF", "ALT"])
    allele_group = melt_df.groupby(["CHROM", "POS", "REF", "ALT"])["value"]
    cal_miss_rate = partial(allele_stats, stat_type=StatType.MISS)
    cal_het_rate = partial(allele_stats, stat_type=StatType.HET)
    cal_maf_rate = partial(allele_stats, stat_type=StatType.MAF)
    miss_rate_df = allele_group.apply(cal_miss_rate)
    het_rate_df = allele_group.apply(cal_het_rate)
    maf_df = allele_group.apply(cal_maf_rate)
    stats_df = pd.concat([miss_rate_df, het_rate_df, maf_df], axis=1)
    stats_df.columns = ["missing", "het", "maf"]
    stats_df.reset_index(inplace=True)
    stats_df["alleles"] = stats_df.apply(lambda x: f'{x["REF"]},{x["ALT"]}', axis=1)
    return stats_df


def vcfStats(gt_table: Path, vcf_stats: Path) -> None:
    dfs = pd.read_table(gt_table, chunksize=100_000)
    for n, df in enumerate(dfs):
        mode = "w" if n == 0 else "a"
        stats_df = transform_one(df)
        stats_df.to_csv(
            vcf_stats,
            mode=mode,
            index=False,
            header=False,
            columns=["CHROM", "POS", "alleles", "missing", "het", "maf"],
        )


if __name__ == "__main__":
    typer.run(vcfStats)
