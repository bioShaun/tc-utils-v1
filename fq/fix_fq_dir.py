#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from enum import StrEnum  # 确保正确导入StrEnum
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import typer
from loguru import logger
from tqdm import tqdm

__version__ = "1.0"

app = typer.Typer(help="FASTQ文件处理和合并工具")

# 配置常量
DEFAULT_BASE_DIR = Path("/public/home/zxchen/data_trans")
FASTQ_EXTENSIONS = ("*.fastq.gz", "*.fq.gz")
READ_TYPE_PATTERNS = {
    "combined_R1.fastq.gz": "R1",
    "combined_R2.fastq.gz": "R2",
    "_R1.fastq.gz": "R1",
    "_R2.fastq.gz": "R2",
    "_1.fastq.gz": "R1",
    "_2.fastq.gz": "R2",
}

# TODO 合并后的表格也需要检查是否有重复


class FastqErrorType(StrEnum):
    """FASTQ文件错误类型"""

    DUPLICATED = "DUPLICATED"
    INCOMPLETE = "INCOMPLETE"
    FORMAT = "FORMAT"
    DATA_SIZE = "DATA_SIZE"
    LOW_DATA = "LOW_DATA"  # 添加缺失的LOW_DATA类型


@dataclass
class FastqError:
    """FASTQ文件错误"""

    name: str
    error_type: FastqErrorType  # 使用字符串类型存储错误类型
    error_message: str


@dataclass
class FastqErrorRecorder:
    """FASTQ文件错误处理器"""

    _errors: List[FastqError] = field(default_factory=list)

    def record_error(self, name: str, error_type: FastqErrorType, error_message: str):
        """记录错误"""
        self._errors.append(FastqError(name, error_type, error_message))

    def get_errors(self) -> List[FastqError]:
        """获取所有错误"""
        return self._errors


class FastqProcessor:
    """FASTQ文件处理器"""

    def __init__(
        self, base_dir: Path, line: str, error_recorder: FastqErrorRecorder
    ):  # 使用Optional类型提示
        self.fq_line_dir = base_dir / line
        self.error_recorder = error_recorder

    def parse_fastq_filename(self, sample_path: Path) -> List[Dict]:
        """解析FASTQ文件名，支持多种命名格式"""
        if not sample_path.exists():
            logger.warning(f"样品路径不存在: {sample_path}")
            return []

        filename = sample_path.name
        fastqs = []

        # 搜索所有可能的FASTQ文件扩展名
        for pattern in FASTQ_EXTENSIONS:
            fastqs.extend(sample_path.glob(pattern))

        lib_id = filename.split("-")[-1]
        if lib_id.isdigit():
            lib_id = "-".join(filename.split("-")[-2:])
        lib_info = []

        if not fastqs:
            logger.warning(f"在 {sample_path} 中未找到FASTQ文件")
            self.error_recorder.record_error(
                lib_id,
                FastqErrorType.INCOMPLETE,  # 使用value属性
                f"未找到FASTQ文件 {sample_path}",
            )
            return []

        for fq in fastqs:
            read_type = self._determine_read_type(fq.name)
            if read_type:
                lib_info.append(
                    {
                        "libid": lib_id,
                        "read_type": read_type,
                        "path": str(fq.absolute()),
                    }
                )
            else:
                logger.warning(f"无法识别的FASTQ文件: {fq.name}")
                self.error_recorder.record_error(
                    name=lib_id,
                    error_type=FastqErrorType.FORMAT,  # 使用value属性
                    error_message=f"无法识别的FASTQ文件: {fq.name}",
                )

        r1_count = sum(1 for fq in fastqs if fq.name.endswith("_R1.fastq.gz"))
        r2_count = sum(1 for fq in fastqs if fq.name.endswith("_R2.fastq.gz"))

        if r1_count != r2_count:
            self.error_recorder.record_error(
                name=lib_id,
                error_type=FastqErrorType.FORMAT,
                error_message=f"R1和R2数量不一致: {r1_count} != {r2_count}",
            )

        return lib_info

    def _libid_not_duplicated(self, fq_line: Path) -> None:
        try:
            sample_dirs = list(fq_line.glob("Sample*"))
        except Exception as e:
            raise Exception(f"遍历目录失败: {e}")

        libids = [dir.name.split("-")[-1] for dir in sample_dirs]
        libid_counts = Counter(libids)
        duplicated = {libid for libid, count in libid_counts.items() if count > 1}

        if duplicated:
            # raise Exception(f"发现重复的 libid: {duplicated}，请检查目录命名！")
            for lib_id in duplicated:
                self.error_recorder.record_error(
                    name=lib_id,
                    error_type=FastqErrorType.DUPLICATED,
                    error_message=f"{fq_line}发现重复的libid: {lib_id}",
                )

    def _determine_read_type(self, filename: str) -> Optional[str]:
        """根据文件名确定读取类型(R1/R2)"""
        for pattern, read_type in READ_TYPE_PATTERNS.items():
            if filename.endswith(pattern):
                return read_type
        return None

    def build_libid_fastq_map(self, fastq_path: Path) -> pd.DataFrame:
        """构建library ID到FASTQ文件的映射"""
        if not fastq_path.exists():
            logger.error(f"FASTQ路径不存在: {fastq_path}")
            return pd.DataFrame()

        libid_map = []
        try:
            sample_dirs = list(fastq_path.glob("Sample*"))
        except Exception as e:
            logger.error(f"遍历目录失败: {e}")
            raise ValueError(f"遍历目录失败: {e}")

        if not sample_dirs:
            logger.warning(f"在 {fastq_path} 中未找到Sample*目录")
            raise ValueError(f"在 {fastq_path} 中未找到Sample*目录")

        for path in tqdm(sample_dirs, desc=f"处理 {fastq_path.name}"):
            try:
                info = self.parse_fastq_filename(path)
                libid_map.extend(info)
            except Exception as e:
                logger.error(f"处理 {path} 时出错: {e}")
                raise ValueError(f"处理 {path} 时出错: {e}")

        return pd.DataFrame(libid_map)

    def build_config(self, force_rebuild: bool = False) -> None:
        """读取或构建配置文件"""
        config_file = self.fq_line_dir / "libid_fastq_config.tsv"

        try:
            logger.info(f"使用已有配置文件: {config_file}")
            df = pd.read_table(config_file)
            # 验证必要的列是否存在
            required_cols = ["libid", "read_type", "path"]
            if not all(col in df.columns for col in required_cols):
                logger.warning(f"配置文件缺少必要列，重新构建: {config_file}")
                raise ValueError("配置文件格式不正确")

        except Exception as e:
            logger.warning(f"读取配置文件失败，重新构建: {e}")


@app.command()
def fix(
    sample_info: Path = typer.Argument(
        ..., help="样品信息TSV文件，必须包含libid、sample_id、dir_name列"
    ),
    base_dir: Path = typer.Option(DEFAULT_BASE_DIR, help="包含所有FASTQ数据的基础目录"),
    output_dir: Optional[Path] = typer.Option(
        Path("raw_data"), help="FASTQ文件输出目录"
    ),
    check_file: Path = typer.Option("check_file.tsv", help="检查结果输出文件"),
    threads: int = typer.Option(8, min=1, max=32, help="并行处理线程数"),
    force_rebuild: bool = typer.Option(False, help="强制重建配置文件"),
    rm_empty_data: bool = typer.Option(True, help="删除空数据文件"),
    empty_data_threshold: int = typer.Option(0.01, help="空数据阈值"),
):
    """
    FASTQ文件处理和合并工具

    此工具用于根据样品信息文件查找、验证和合并FASTQ文件。
    """

    # 验证样品信息
    # 初始化错误收集器
    error_collector = FastqErrorRecorder()
    processor = FastqProcessor(base_dir, error_recorder=error_collector)


if __name__ == "__main__":
    app()
