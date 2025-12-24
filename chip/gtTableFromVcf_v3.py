import re
from enum import Enum
from functools import partial
from pathlib import Path

import delegator
import pandas as pd
import typer
from loguru import logger
from pandarallel import pandarallel

LOCATION_COLS = ["CHROM", "POS", "REF", "ALT"]


class OutType(str, Enum):
    txt = "txt"
    xlsx = "xlsx"


class TableColumn(Enum):
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


class GT_VALUE(Enum):
    NA = "./."
    HET = "0/1"
    REF = "0/0"
    ALT = "1/1"


class MissFmt(str, Enum):
    N = "N"
    NN = "NN"
    DASH = "--"

    def __str__(self) -> str:
        return self.value


def vcf2gt(
    vcf_file: Path, target_id: Path, out_file: Path, threads: int, force: bool = False
) -> Path:
    gt_file = Path(f"{out_file}.gt.txt.gz")
    if gt_file.exists() and not force:
        return gt_file
    gt_file.parent.mkdir(parents=True, exist_ok=True)
    add_id_vcf = vcf_file.parent / f"add-id.{vcf_file.name}"
    target_vcf = f"{out_file}.vcf.gz"
    logger.info("add id to vcf")
    cmd1 = f'bcftools annotate --set-id "%CHROM\\_%POS" {vcf_file} -Oz -o {add_id_vcf} --threads {threads}'
    delegator.run(cmd1)
    logger.info("extract target vcf")
    cmd2 = f"bcftools view -i ID=@{target_id} {add_id_vcf} -Oz -o {target_vcf} --threads {threads}"
    logger.info(f"run: {cmd2}")
    delegator.run(cmd2)
    logger.info("extract target genotype")
    cmd3 = f'bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n" {target_vcf} | sed -re "s;\\|;/;g" | gzip > {gt_file}'
    logger.info(f"run: {cmd3}")
    delegator.run(cmd3)
    return gt_file


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
    # melt_gt_df[TableColumn.GENOTYPE.value] = melt_gt_df.apply(myConvertGT, axis=1)

    convert_df = melt_gt_df.set_index(
        [*LOCATION_COLS, TableColumn.SAMPLE_NAME.value]
    ).unstack(4)
    convert_df.columns = convert_df.columns.droplevel()

    out_cols = [each for each in gt_df.columns if each not in [*LOCATION_COLS]]
    out_df = convert_df[out_cols].reset_index()
    return loc_df.merge(out_df)


def transformOneAlt(ref: str, alt: str) -> str:
    """
    将VCF变异转换为注释格式
    - SNP: 只返回变化的碱基
    - Insertion/Deletion: 返回插入或删除的序列
    - 遇到delins复杂变异: 抛出异常
    """
    if not ref or not alt:
        raise ValueError("ref and alt cannot be empty")
    
    if alt == "*":
        return f"del{ref}"
    
    # 去除共同前缀
    prefix_len = 0
    min_len = min(len(ref), len(alt))
    while prefix_len < min_len and ref[prefix_len] == alt[prefix_len]:
        prefix_len += 1
    
    # 去除共同后缀
    suffix_len = 0
    while (suffix_len < min_len - prefix_len and 
           ref[-(suffix_len + 1)] == alt[-(suffix_len + 1)]):
        suffix_len += 1
    
    # 提取变化部分
    if suffix_len > 0:
        ref_changed = ref[prefix_len:-suffix_len]
        alt_changed = alt[prefix_len:-suffix_len]
    else:
        ref_changed = ref[prefix_len:]
        alt_changed = alt[prefix_len:]
    
    # 纯删除
    if ref_changed and not alt_changed:
        return f"del{ref_changed}"
    
    # 纯插入
    if not ref_changed and alt_changed:
        return f"ins{alt_changed}"
    
    # SNP (单碱基替换)
    if len(ref_changed) == 1 and len(alt_changed) == 1:
        return alt_changed
    
    # MNP (多碱基替换,长度相等)
    if len(ref_changed) == len(alt_changed):
        return alt_changed
    
    # 复杂变异 (delins) - 报错
    raise ValueError(
        f"Complex variant (delins) detected: "
        f"ref={ref}, alt={alt}, "
        f"ref_changed={ref_changed}, alt_changed={alt_changed}"
    )


def transformAlt(ref: str, alts: str):
    alt_list = alts.split(",")
    new_alt_list = [transformOneAlt(ref, each) for each in alt_list]
    return ",".join(new_alt_list)


def get_sample_names(vcf_file: Path) -> list:
    cmd = f"bcftools query -l {vcf_file}"
    logger.info(f"run: {cmd}")
    return delegator.run(cmd).out.strip().split("\n")


def main(
    vcf: Path,
    target_id: Path,
    out_file: Path,
    miss_fmt: str = "NN",
    threads: int = 4,
    gt_sep: str = "",
):
    pandarallel.initialize(progress_bar=True, nb_workers=threads)
    gt_file = vcf2gt(vcf, target_id, out_file, threads)
    sample_list = get_sample_names(vcf)
    columns = LOCATION_COLS + sample_list
    gt_dfs = pd.read_csv(
        gt_file,
        sep="\t",
        header=None,
        names=columns,
        chunksize=10000,
    )
    # va_type_df = pd.read_table(va_type)
    for i, gt_df in enumerate(gt_dfs):
        gt_df["ALT"] = gt_df.parallel_apply(
            lambda x: transformAlt(x.REF, x.ALT), axis=1
        )
        gt_df["REF"] = gt_df["REF"].map(lambda x: x[0])
        gt_df = gt_df.set_index(LOCATION_COLS)
        gt_df.replace(".", "./.", inplace=True)
        gt_df.replace("1/0", "0/1", inplace=True)
        gt_df = gt_df.reset_index()
        logger.info(f"Transforming {i}")
        seq_df = gt2seq(gt_df, miss_fmt, gt_sep)
        # gt_df = va_type_df.merge(gt_df)
        # seq_df = va_type_df.merge(seq_df)
        mode = "w"
        header = True
        if i > 0:
            mode = "a"
            header = False
        gt_df.to_csv(
            f"{out_file}.gt.txt.gz", index=False, sep="\t", header=header, mode=mode
        )
        seq_df.to_csv(
            f"{out_file}.seq.txt.gz", sep="\t", header=header, mode=mode, index=False
        )


if __name__ == "__main__":
    typer.run(main)
