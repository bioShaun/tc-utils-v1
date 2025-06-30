from enum import StrEnum
from functools import partial
from pathlib import Path
from typing import Dict, List

import delegator
import pandas as pd
import typer
from loguru import logger
from pandarallel import pandarallel
from tqdm import tqdm
from typing_extensions import Tuple

LOCATION_COLS = ["CHROM", "POS", "REF", "ALT"]


class GT_VALUE(StrEnum):
    NA = "./."
    HET = "0/1"
    REF = "0/0"
    ALT = "1/1"


class VA_TYPE(StrEnum):
    SNP = "SNP"
    INS = "INS"
    DEL = "DEL"
    MIXED = "MIXED"


class TableColumn(StrEnum):
    CHROM = "CHROM"
    POS = "POS"
    REF = "REF"
    ALT = "ALT"
    SAMPLE_NAME = "SAMPLE_NAME"
    GENOTYPE = "GENOTYPE"
    P1 = "P1"
    P2 = "P2"
    STATS = "stats"
    COUNT = "Count"
    RATIO = "Ratio"


def vcf2gt(vcf_file: Path, force: bool = False) -> Path:
    gt_file = vcf_file.with_suffix(".gt.txt.gz")
    if gt_file.exists() and not force:
        return gt_file
    gt_file.parent.mkdir(parents=True, exist_ok=True)
    cmd = f'bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n" {vcf_file} | sed -re "s;\\|;/;g" | gzip > {gt_file}'
    logger.info(f"run: {cmd}")
    delegator.run(cmd)
    return gt_file


def get_sample_names(vcf_file: Path) -> list:
    cmd = f"bcftools query -l {vcf_file}"
    logger.info(f"run: {cmd}")
    return delegator.run(cmd).out.strip().split("\n")


def norm_alleles(vcf_file: Path, ref_fa: Path, force: bool = False) -> Path:
    norm_file = vcf_file.with_suffix(".norm.vcf.gz")
    if norm_file.exists() and not force:
        return norm_file
    cmd = f'bcftools norm -f {ref_fa} -m -any  {vcf_file} | bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\n" | gzip > {norm_file}'
    logger.info(f"run: {cmd}")
    delegator.run(cmd)
    return norm_file


def transformOneAlt(ref: str, alt: str) -> str:
    if alt == "*":
        return f"del{ref}"
    if len(ref) > len(alt):
        del_part = ref[1:]
        return f"del{del_part}"
    if len(ref) < len(alt):
        ins_part = alt[1:]
        return f"ins{ins_part}"
    return alt


def load_allele_map_dict(norm_gt: pd.DataFrame) -> Dict[Tuple, List[str]]:
    out_dict = {}
    for row in tqdm(norm_gt.itertuples(), total=len(norm_gt)):
        va_id = (row.CHROM, row.POS)
        out_dict.setdefault(va_id, []).append(transformOneAlt(row.REF, row.ALT))
    return out_dict


def npConvertGT(row: pd.Series, miss_fmt: str, gt_sep: str) -> str:
    if row[TableColumn.GENOTYPE.value] == GT_VALUE.NA.value:
        return miss_fmt
    allele1, allele2 = [each for each in row[TableColumn.GENOTYPE.value].split("/")]
    allele_list = [
        row[TableColumn.REF.value],
        *row[TableColumn.ALT.value].split(","),
    ]
    allele1_seq = "N" if allele1 == "." else allele_list[int(allele1)]
    allele2_seq = "N" if allele2 == "." else allele_list[int(allele2)]
    if allele1 == allele2:
        if len(allele1_seq) > 1:
            return allele1_seq
        return f"{allele1_seq}{gt_sep}{allele2_seq}"
    if len(allele1_seq) > 1 or len(allele2_seq) > 1:
        return f"{allele1_seq}/{allele2_seq}"
    return f"{allele1_seq}{gt_sep}{allele2_seq}"


def gt2seq(gt_df: pd.DataFrame, miss_fmt: str, gt_sep: str):
    # transform
    gt_df.drop_duplicates(subset=LOCATION_COLS, inplace=True)
    loc_df = gt_df[LOCATION_COLS].copy()
    melt_gt_df = gt_df.melt(
        id_vars=LOCATION_COLS,
        var_name=TableColumn.SAMPLE_NAME.value,
        value_name=TableColumn.GENOTYPE.value,
    )
    myConvertGT = partial(npConvertGT, miss_fmt=miss_fmt, gt_sep=gt_sep)
    melt_gt_df[TableColumn.GENOTYPE.value] = melt_gt_df.parallel_apply(
        myConvertGT, axis=1
    )  # type: ignore

    convert_df = melt_gt_df.set_index(
        [*LOCATION_COLS, TableColumn.SAMPLE_NAME.value]
    ).unstack(4)
    convert_df.columns = convert_df.columns.droplevel()

    out_cols = [each for each in gt_df.columns if each not in [*LOCATION_COLS]]
    out_df = convert_df[out_cols].reset_index()
    return loc_df.merge(out_df)


def variant_type(alt: str) -> str:
    type_list = []
    for each_alt in alt.split(","):
        if "ins" in each_alt:
            type_list.append(VA_TYPE.INS.value)
        elif "del" in each_alt:
            type_list.append(VA_TYPE.DEL.value)
        else:
            type_list.append(VA_TYPE.SNP.value)
    if len(set(type_list)) == 1:
        return type_list[0]
    else:
        return VA_TYPE.MIXED.value


def main(
    vcf: Path,
    ref_fa: Path,
    out_prefix: Path,
    miss_fmt: str = "NN",
    threads: int = 4,
    gt_sep: str = "",
) -> None:
    pandarallel.initialize(progress_bar=True, nb_workers=threads)
    gt_file = vcf2gt(vcf)
    sample_list = get_sample_names(vcf)
    gt_dfs = pd.read_table(
        gt_file,
        header=None,
        names=[*LOCATION_COLS, *sample_list],
        chunksize=10000,
    )

    norm_allele_file = norm_alleles(vcf, ref_fa)
    logger.info(f"left align alleles")
    norm_gt_df = pd.read_table(
        norm_allele_file,
        header=None,
        names=[*LOCATION_COLS],
    )
    allele_dict = load_allele_map_dict(norm_gt_df)

    for i, gt_df in enumerate(gt_dfs):
        gt_df["REF"] = gt_df["REF"].map(lambda x: x[0])
        gt_df["ALT"] = gt_df.apply(
            lambda x: ",".join(allele_dict[(x.CHROM, x.POS)]), axis=1
        )
        gt_df = gt_df.set_index(LOCATION_COLS)
        gt_df.replace(".", "./.", inplace=True)
        gt_df.replace("1/0", "0/1", inplace=True)
        gt_df = gt_df.reset_index()
        logger.info(f"Transforming {i}")
        seq_df = gt2seq(gt_df, miss_fmt, gt_sep)
        gt_df["TYPE"] = gt_df["ALT"].map(variant_type)
        seq_df["TYPE"] = seq_df["ALT"].map(variant_type)
        mode = "w"
        header = True
        if i > 0:
            mode = "a"
            header = False
        gt_df.to_csv(
            f"{out_prefix}.gt.txt.gz",
            index=False,
            sep="\t",
            header=header,
            mode=mode,
            columns=[*LOCATION_COLS, "TYPE", *sample_list],
        )
        seq_df.to_csv(
            f"{out_prefix}.seq.txt.gz",
            sep="\t",
            header=header,
            mode=mode,
            index=False,
            columns=[*LOCATION_COLS, "TYPE", *sample_list],
        )


if __name__ == "__main__":
    typer.run(main)
