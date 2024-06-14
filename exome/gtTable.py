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
    return f"{allele_list[allele1]}{allele_list[allele2]}".replace("*", "N")


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


def main(
    gt_file: Path,
    sample_file: Path,
    out_file: Path,
    miss_fmt: str = "NN",
    threads: int = 4,
    xlsx: bool = False,
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
    seq_dfs = []
    gt_concat_dfs = []
    for i, gt_df in enumerate(gt_dfs):
        gt_df = gt_df.set_index(LOCATION_COLS)
        gt_df.replace(".", "./.", inplace=True)
        gt_df.replace("1/0", "0/1", inplace=True)
        gt_df = gt_df.reset_index()
        logger.info(f"Transforming {i}")
        seq_df = gt2seq(gt_df, miss_fmt)

        # seq_df = pd.concat(seq_dfs)
        # gt_df = pd.concat(gt_concat_dfs)
        mode = "w"
        header = True
        if i > 0:
            mode = "a"
            header = False
        if xlsx:
            seq_dfs.append(seq_df)
            gt_concat_dfs.append(gt_df)
        else:
            gt_df.to_csv(
                f"{out_file}.gt.txt.gz", index=False, sep="\t", header=header, mode=mode
            )
            seq_df.to_csv(
                f"{out_file}.seq.txt.gz",
                sep="\t",
                header=header,
                mode=mode,
                index=False,
            )

    if xlsx:
        merged_gt_df = pd.concat(gt_concat_dfs)
        merged_seq_df = pd.concat(seq_dfs)
        merged_gt_df.to_excel(f"{out_file}.gt.xlsx", index=False)
        merged_seq_df.to_excel(f"{out_file}.seq.xlsx", index=False)


if __name__ == "__main__":
    typer.run(main)
