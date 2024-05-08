import re
from enum import Enum
from functools import partial
from pathlib import Path

import pandas as pd
import typer
from loguru import logger
from pandarallel import pandarallel

LOCATION_COLS = ["CHROM", "POS", "REF", "ALT"]


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


def npConvertGT(row: pd.Series, miss_fmt: str) -> str:
    if row[TableColumn.GENOTYPE.value] == GT_VALUE.NA.value:
        return miss_fmt
    allele1, allele2 = [
        int(each) for each in row[TableColumn.GENOTYPE.value].split("/")
    ]
    allele_list = [
        row[TableColumn.REF.value],
        *row[TableColumn.ALT.value].split(","),
    ]
    if allele1 == allele2:
        if len(allele_list[allele1]) > 1:
            return allele_list[allele1]
        return f"{allele_list[allele1]}{allele_list[allele2]}"
    if len(allele_list[allele1]) > 1 or len(allele_list[allele2]) > 1:
        return f"{allele_list[allele1]}/{allele_list[allele2]}"
    return f"{allele_list[allele1]}{allele_list[allele2]}"


def gt2seq(gt_df: pd.DataFrame, miss_fmt: str):
    # transform
    gt_df.drop_duplicates(subset=LOCATION_COLS, inplace=True)
    loc_df = gt_df[LOCATION_COLS].copy()
    melt_gt_df = gt_df.melt(
        id_vars=LOCATION_COLS,
        var_name=TableColumn.SAMPLE_NAME.value,
        value_name=TableColumn.GENOTYPE.value,
    )
    myConvertGT = partial(npConvertGT, miss_fmt=miss_fmt)
    melt_gt_df[TableColumn.GENOTYPE.value] = melt_gt_df.parallel_apply(
        myConvertGT, axis=1
    )
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


def main(
    gt_file: Path,
    sample_file: Path,
    # va_type: Path,
    out_file: Path,
    miss_fmt: str = "NN",
    threads: int = 4,
):
    pandarallel.initialize(progress_bar=True, nb_workers=threads)
    sample_list = pd.read_csv(sample_file, header=None)[0].tolist()
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
        seq_df = gt2seq(gt_df, miss_fmt)
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
