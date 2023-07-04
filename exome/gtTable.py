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
    GT_MAP = {
        GT_VALUE.ALT.value: f"{row[TableColumn.ALT.value]}{row[TableColumn.ALT.value]}",
        GT_VALUE.REF.value: f"{row[TableColumn.REF.value]}{row[TableColumn.REF.value]}",
        GT_VALUE.HET.value: f"{row[TableColumn.REF.value]}{row[TableColumn.ALT.value]}",
        GT_VALUE.NA.value: miss_fmt,
    }
    return GT_MAP[row[TableColumn.GENOTYPE.value]]


def gt2seq(gt_df: pd.DataFrame, miss_fmt: str):
    # transform
    melt_gt_df = gt_df.melt(
        id_vars=LOCATION_COLS,
        var_name=TableColumn.SAMPLE_NAME.value,
        value_name=TableColumn.GENOTYPE.value,
    )
    myConvertGT = partial(npConvertGT, miss_fmt=miss_fmt)
    melt_gt_df[TableColumn.GENOTYPE.value] = melt_gt_df.parallel_apply(
        myConvertGT, axis=1
    )

    convert_df = melt_gt_df.set_index(
        [*LOCATION_COLS, TableColumn.SAMPLE_NAME.value]
    ).unstack(4)
    convert_df.columns = convert_df.columns.droplevel()

    out_cols = [each for each in gt_df.columns if each not in [*LOCATION_COLS]]
    return convert_df[out_cols]


def main(
    gt_file: Path,
    sample_file: Path,
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
    seq_dfs = []
    gt_concat_dfs = []
    for i, gt_df in enumerate(gt_dfs):
        logger.info(f"Transforming {i}")
        seq_df = gt2seq(gt_df, miss_fmt)
        seq_dfs.append(seq_df)
        gt_concat_dfs.append(gt_df)
    seq_df = pd.concat(seq_dfs)
    gt_df = pd.concat(gt_concat_dfs)
    gt_df.to_csv(f"{out_file}.gt.txt.gz", index=False, sep="\t")
    seq_df.to_csv(f"{out_file}.seq.txt.gz", sep="\t")


if __name__ == "__main__":
    typer.run(main)
