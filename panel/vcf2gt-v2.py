from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import List, Optional

import delegator
import numpy as np
import pandas as pd
import typer
from loguru import logger

LOCATION_COLS = ["CHROM", "POS", "REF", "ALT"]


@dataclass
class Config:
    CHUNK_SIZE: int = 10000
    DEFAULT_MISS_FMT: str = "NN"
    DEFAULT_GT_SEP: str = ""


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


# exceptions.py
class VCFProcessError(Exception):
    """VCF处理相关的自定义异常"""

    pass


class DataValidationError(VCFProcessError):
    """数据验证错误"""

    pass


class FileProcessError(VCFProcessError):
    """文件处理错误"""

    pass


def validate_input_file(file_path: Path) -> bool:
    """验证输入文件"""
    if not file_path.exists():
        raise FileNotFoundError(f"Input file not found: {file_path}")
    return True


def vcf2gt(vcf_file: Path) -> Path:
    gt_file = vcf_file.with_suffix(".gt.txt.gz")
    gt_file.parent.mkdir(parents=True, exist_ok=True)
    cmd = f'bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n" {vcf_file} | sed -re "s;\\|;/;g" | gzip > {gt_file}'
    logger.info(f"run: {cmd}")
    delegator.run(cmd)
    return gt_file


def get_sample_names(vcf_file: Path) -> list:
    cmd = f"bcftools query -l {vcf_file}"
    logger.info(f"run: {cmd}")
    return delegator.run(cmd).out.strip().split("\n")


class VCFProcessor:
    """VCF文件处理类"""

    def __init__(
        self,
        input_file: Path,
        output_prefix: Path,
        input_type: InputType = InputType.VCF,
        miss_fmt: str = Config.DEFAULT_MISS_FMT,
        gt_sep: str = Config.DEFAULT_GT_SEP,
        sample_file: Optional[Path] = None,
    ):
        self.input_file = input_file
        self.output_prefix = output_prefix
        self.input_type = input_type
        self.miss_fmt = miss_fmt
        self.gt_sep = gt_sep
        self.sample_file = sample_file
        self.gt_df = None

    def load_data(self) -> None:
        """加载数据"""
        logger.info(f"Loading data from {self.input_file}")
        try:
            if self.input_type == InputType.VCF:
                self.gt_df = self._load_vcf()
            else:
                self.gt_df = self._load_table()
        except Exception as e:
            raise DataValidationError(f"Failed to load data: {str(e)}")

    def _load_vcf(self) -> pd.DataFrame:
        """加载VCF文件"""
        gt_file = vcf2gt(self.input_file)
        sample_list = get_sample_names(self.input_file)
        return pd.read_table(gt_file, header=None, names=[*LOCATION_COLS, *sample_list])

    def _load_table(self) -> pd.DataFrame:
        """加载表格文件"""
        if self.sample_file is None:
            raise DataValidationError("Sample file is required for TABLE input type")
        sample_list = pd.read_csv(self.sample_file, header=None)[0].tolist()
        return pd.read_table(
            self.input_file, header=None, names=[*LOCATION_COLS, *sample_list]
        )

    def process_data(self) -> None:
        """处理数据"""
        logger.info("Processing data")
        if self.gt_df is None:
            raise DataValidationError("No data loaded")
        self.gt_df = self.gt2seq(self.gt_df)

    def gt2seq(self, gt_df: pd.DataFrame) -> pd.DataFrame:
        """基因型转换为序列"""
        try:
            # 数据预处理
            gt_df = gt_df.copy()
            gt_df.drop_duplicates(subset=LOCATION_COLS, inplace=True)

            # 向量化转换
            def vectorized_convert_gt(genotype, ref, alt):
                mask = genotype == GT_VALUE.NA.value
                allele_list = np.array([ref] + alt.split(","))
                allele1, allele2 = np.char.split(genotype, "/")
                return np.where(
                    mask,
                    self.miss_fmt,
                    np.where(
                        allele1 == allele2,
                        allele_list[allele1.astype(int)],
                        allele_list[allele1.astype(int)]
                        + self.gt_sep
                        + allele_list[allele2.astype(int)],
                    ),
                )

            # 应用转换
            for col in gt_df.columns:
                if col not in LOCATION_COLS:
                    gt_df[col] = vectorized_convert_gt(
                        gt_df[col], gt_df[TableColumn.REF], gt_df[TableColumn.ALT]
                    )

            return gt_df

        except Exception as e:
            raise VCFProcessError(f"Failed to convert genotypes: {str(e)}") from e

    def save_results(self) -> None:
        """保存结果"""
        logger.info(f"Saving results to {self.output_prefix}")
        try:
            output_file = self.output_prefix.with_suffix(".csv")
            self.gt_df.to_csv(output_file, index=False)
        except Exception as e:
            raise FileProcessError(f"Failed to save results: {str(e)}") from e


# main.py


def main(
    input_file: Path,
    out_prefix: Path,
    input_type: InputType = InputType.VCF,
    force: bool = False,
    sample_file: Optional[Path] = typer.Option(None),
    miss_fmt: str = "NN",
    gt_sep: str = "",
) -> None:
    """
    处理VCF文件或基因型表格文件
    """
    try:
        # 设置日志
        logger.info(f"Starting processing file: {input_file}")

        # 验证输入
        validate_input_file(input_file)
        if input_type == InputType.TABLE and sample_file:
            validate_input_file(sample_file)

        # 处理数据
        processor = VCFProcessor(
            input_file=input_file,
            output_prefix=out_prefix,
            input_type=input_type,
            miss_fmt=miss_fmt,
            gt_sep=gt_sep,
            sample_file=sample_file,
        )

        processor.load_data()
        processor.process_data()
        processor.save_results()

        logger.success("Processing completed successfully")

    except Exception as e:
        logger.error(f"Processing failed: {str(e)}")
        raise VCFProcessError(f"Failed to process file: {str(e)}") from e


if __name__ == "__main__":
    typer.run(main)
