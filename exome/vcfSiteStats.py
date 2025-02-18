import re
from pathlib import Path
from typing import Tuple

import delegator
import pandas as pd
import typer
from pandarallel import pandarallel
from tqdm import tqdm


def get_gt_type(gt: str) -> str:
    if gt in ["./.", "."]:
        return "miss"
    try:
        allele1, allele2 = gt.split("/")[:2]
    except ValueError:
        print("error gt:", gt)
        raise ValueError(f"{gt} format error!")
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
    df["indel_length"] = df.parallel_apply(get_index_len, axis=1)  # type: ignore
    df["ALT"] = df.parallel_apply(lambda x: transformAlt(x.REF, x.ALT), axis=1)  # type: ignore
    df["REF"] = df["REF"].map(lambda x: x[0])
    df["id"] = df.apply(lambda x: f"{x['CHROM']}_{x['POS']}", axis=1)
    df["stats_info"] = df.parallel_apply(allele_stats, axis=1)  # type: ignore
    df[["alleles", "missing", "het", "maf"]] = pd.DataFrame(
        df["stats_info"].tolist(), index=df.index
    )
    df["indel_type"] = df["alleles"].map(get_indel_type)
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
    ref_len = len(row.REF)
    alt_list = [len(x) for x in row.ALT.split(",")]
    alt_len_max, alt_len_min = max(alt_list), min(alt_list)
    return max(abs(ref_len - alt_len_max), abs(ref_len - alt_len_min))


def get_indel_type(alleles: str) -> str:
    if "ins" in alleles:
        if "del" in alleles:
            return "MIXED"
        return "INS"
    elif "del" in alleles:
        return "DEL"
    return "SNP"


def vcf2gt(vcf_file: Path, force: bool = False) -> Path:
    gt_file = vcf_file.with_suffix(".gt.txt.gz")
    if gt_file.exists() and not force:
        return gt_file
    gt_file.parent.mkdir(parents=True, exist_ok=True)
    cmd = f'bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n" {vcf_file} | sed -re "s;\\|;/;g" | gzip > {gt_file}'
    delegator.run(cmd)
    return gt_file


def indel_right_pos(row: pd.Series) -> int:
    if row["indel_type"] in ["DEL", "MIXED"]:
        return row["pos"] + row["indel_length"]
    return row["pos"]


def vcfStats(
    vcf: Path, vcf_stats: Path, threads: int = 4, force: bool = False, snp: bool = True
) -> None:
    pandarallel.initialize(nb_workers=threads)
    gt_table = vcf2gt(vcf_file=vcf, force=force)
    dfs = pd.read_table(gt_table, chunksize=100_000, header=None)
    left_file = vcf_stats.with_suffix(".left.tsv")
    right_file = vcf_stats.with_suffix(".right.tsv")
    out_cols = [
        "id",
        "chrom",
        "pos",
        "alleles",
        "missing",
        "het",
        "maf",
        "indel_length",
        "indel_type",
    ]
    if snp:
        out_cols = ["id", "chrom", "pos", "alleles", "missing", "het", "maf"]
    for n, df in tqdm(enumerate(dfs)):
        mode = "w" if n == 0 else "a"
        header = True if n == 0 else False
        stats_df = transform_one(df)
        out_df = stats_df[
            [
                "id",
                "CHROM",
                "POS",
                "alleles",
                "missing",
                "het",
                "maf",
                "indel_length",
                "indel_type",
            ]
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
            "indel_type",
        ]
        # choose best maf one for multi alt sites
        out_df = out_df.loc[out_df.groupby(["chrom", "pos"])["maf"].idxmax()]
        out_df.to_csv(
            vcf_stats,
            mode=mode,
            index=False,
            float_format="%.3f",
            sep="\t",
            header=header,
            columns=out_cols,
        )
        if not snp:
            out_df.to_csv(
                left_file,
                mode=mode,
                index=False,
                float_format="%.3f",
                sep="\t",
                header=header,
            )
            right_df = out_df.copy()
            right_df["pos"] = right_df.apply(indel_right_pos, axis=1)
            right_df.to_csv(
                right_file,
                mode=mode,
                index=False,
                float_format="%.3f",
                sep="\t",
                header=header,
            )


if __name__ == "__main__":
    typer.run(vcfStats)


