from pathlib import Path
from typing import Tuple

import pandas as pd
import typer
from tqdm import tqdm
from pandarallel import pandarallel


def allele_stats(row: pd.Series) -> Tuple[str, float, float, float]:
    allele_count = row["genotype"].value_counts()
    miss_count = allele_count.get("./.", 0)
    het_count = allele_count.get("0/1", 0)
    real_sample_count = len(row) - miss_count
    miss_rate = miss_count / len(row)
    het_rate = het_count / real_sample_count
    alt_count = allele_count.get("1/1", 0)
    alt_rate = (alt_count * 2 + het_count) / (real_sample_count * 2)
    alleles = f'{row["REF"]},{row["ALT"]}'
    return alleles, miss_rate, het_rate, alt_rate


def transform_one(df: pd.DataFrame) -> pd.DataFrame:
    vcf_columns = ["CHROM", "POS", "REF", "ALT"]
    sample_count = len(df.columns) - len(vcf_columns)
    sample_columns = [f"genotype" for i in range(sample_count)]
    df.columns = [*vcf_columns, *sample_columns]
    df["stats_info"] = df.parallel_apply(allele_stats, axis=1)
    df[["alleles", "missing", "het", "af"]] = pd.DataFrame(
        df["stats_info"].tolist(), index=df.index
    )
    return df


def vcfStats(gt_table: Path, vcf_stats: Path, threads: int = 4) -> None:
    pandarallel.initialize(nb_workers=threads)
    dfs = pd.read_table(gt_table, chunksize=100_000, header=None)
    for n, df in tqdm(enumerate(dfs)):
        mode = "w" if n == 0 else "a"
        stats_df = transform_one(df)
        stats_df.to_csv(
            vcf_stats,
            mode=mode,
            index=False,
            header=False,
            columns=["CHROM", "POS", "alleles", "missing", "het", "af"],
            float_format="%.3f",
        )


if __name__ == "__main__":
    typer.run(vcfStats)
