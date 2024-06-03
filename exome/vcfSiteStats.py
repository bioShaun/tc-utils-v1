import re
from pathlib import Path
from typing import Tuple

import pandas as pd
import typer
from pandarallel import pandarallel
from tqdm import tqdm


def get_gt_type(gt: str) -> str:
    if gt == "./.":
        return "miss"
    allele1, allele2 = gt.split("/")[:2]
    if allele1 != allele2:
        return "het"
    if allele1 == "0":
        return "hom"
    return "alt"


def allele_stats(row: pd.Series) -> Tuple[str, float, float, float]:
    gt_type = row["genotype"].map(get_gt_type)
    allele_count = gt_type.value_counts()
    miss_count = allele_count.get("miss", 0)
    het_count = allele_count.get("het", 0)
    real_sample_count = len(row) - miss_count
    miss_rate = miss_count / len(row)
    het_rate = het_count / real_sample_count
    alt_count = allele_count.get("alt", 0)
    alt_rate = (alt_count * 2 + het_count) / (real_sample_count * 2)
    maf = min(alt_rate, 1 - alt_rate)
    alleles = f'{row["REF"]}/{row["ALT"]}'
    return alleles, miss_rate, het_rate, maf


def transform_one(df: pd.DataFrame) -> pd.DataFrame:
    vcf_columns = ["CHROM", "POS", "REF", "ALT"]
    sample_count = len(df.columns) - len(vcf_columns)
    sample_columns = [f"genotype" for i in range(sample_count)]
    df.columns = [*vcf_columns, *sample_columns]
    df["stats_info"] = df.parallel_apply(allele_stats, axis=1)  # type: ignore
    df[["alleles", "missing", "het", "af"]] = pd.DataFrame(
        df["stats_info"].tolist(), index=df.index
    )
    return df


def transformOneAlt(ref: str, alt: str) -> str:
    if alt == "*":
        return f"del{ref}"
    elif len(ref) > len(alt):
        del_part = re.sub(f"^{alt}", "", ref)
        return f"del{del_part}"
    elif len(ref) < len(alt):
        ins_part = re.sub(f"^{ref}", "", alt)
        return f"ins{ins_part}"
    else:
        if len(ref) == 1:
            return alt
        return alt[0]


def transformAlt(ref: str, alts: str):
    alt_list = alts.split(",")
    new_alt_list = [transformOneAlt(ref, each) for each in alt_list]
    return ",".join(new_alt_list)


def get_index_len(row: pd.Series) -> int:
    ref_len, alt_len = len(row.REF), max(len(x) for x in row.ALT.split(","))
    return abs(ref_len - alt_len)


def vcfStats(gt_table: Path, vcf_stats: Path, threads: int = 4) -> None:
    pandarallel.initialize(nb_workers=threads)
    dfs = pd.read_table(gt_table, chunksize=100_000, header=None)
    for n, df in tqdm(enumerate(dfs)):
        mode = "w" if n == 0 else "a"
        df["indel_length"] = df.parallel_apply(get_index_len, axis=1)
        df["ALT"] = df.parallel_apply(lambda x: transformAlt(x.REF, x.ALT), axis=1)
        df["id"] = df.apply(lambda x: f"{x['CHROM']}_{x['POS']}", axis=1)
        stats_df = transform_one(df)
        out_df = stats_df[
            ["id", "CHROM", "POS", "alleles", "missing", "het", "maf" "indel_length"]
        ].copy()
        out_df.columns = [
            "id",
            "chrom",
            "pos",
            "alleles",
            "missing",
            "het",
            "maf",
            "indel_length",
        ]
        out_df.to_csv(
            vcf_stats,
            mode=mode,
            index=False,
            float_format="%.3f",
        )


if __name__ == "__main__":
    typer.run(vcfStats)
