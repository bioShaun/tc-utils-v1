import re
from enum import StrEnum
from functools import partial
from pathlib import Path

import delegator
import pandas as pd
import typer
from loguru import logger
from pandarallel import pandarallel

LOCATION_COLS = ["CHROM", "POS", "REF", "ALT"]


class InputType(StrEnum):
    VCF = "vcf"
    TABLE = "table"


class TableColumn(StrEnum):
    CHROM = "CHROM"
    POS = "POS"
    REF = "REF"
    ALT = "ALT"
    SAMPLE_NAME = "SAMPLE_NAME"
    GENOTYPE = "GENOTYPE"


class GT_VALUE(StrEnum):
    NA = "./."
    HET = "0/1"
    REF = "0/0"
    ALT = "1/1"


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
    gt_df = gt_df.copy()
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
    if alt == "*":
        return f"del{ref}"
    if len(ref) > len(alt):
        del_part = re.sub(f"^{alt}", "", ref)
        return f"del{del_part}"
    if len(ref) < len(alt):
        ins_part = re.sub(f"^{ref}", "", alt)
        return f"ins{ins_part}"
    if len(ref) == 1:
        return alt
    return alt[0]


def main(
    input_file: Path,
    out_prefix: Path,
    input_type: InputType = InputType.VCF,
    force: bool = False,
    sample_file: Path = typer.Option(None),
    transform_alt: bool = False,
    miss_fmt: str = "NN",
    gt_sep: str = "",
    annotation: Path = typer.Option(None),
) -> None:
    if input_type == InputType.VCF:
        gt_file = vcf2gt(input_file, force=force)
        sample_list = get_sample_names(input_file)
    else:
        gt_file = input_file
        if sample_file is None:
            raise ValueError("Must provide sample_list if input_type is not VCF")
        sample_list = pd.read_csv(sample_file, header=None)[0].tolist()
    gt_df = pd.read_table(gt_file, header=None, names=[*LOCATION_COLS, *sample_list])
    if transform_alt:
        gt_df["ALT"] = gt_df.parallel_apply(
            lambda x: transformAlt(x.REF, x.ALT), axis=1
        )
        gt_df["REF"] = gt_df["REF"].map(lambda x: x[0])
    gt_df = gt_df.set_index(LOCATION_COLS)
    gt_df.replace(".", "./.", inplace=True)
    gt_df.replace("1/0", "0/1", inplace=True)
    gt_df = gt_df.reset_index()
    seq_df = gt2seq(gt_df, miss_fmt, gt_sep)
    seq_df = seq_df.reset_index()
    if annotation:
        anno_df = pd.read_table(annotation)
        gt_df = gt_df.merge(anno_df, how="left")
        seq_df = seq_df.merge(anno_df, how="left")
    gt_df.to_csv(f"{out_prefix}.genotype.01.xlsx", index=False, na_rep="--")
    seq_df.to_csv(f"{out_prefix}.genotype.seq.xlsx", sep="\t", na_rep="--")


if __name__ == "__main__":
    typer.run(main)
