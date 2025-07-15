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

    def record_error(self, name: str, error_type: str, error_message: str):
        """记录错误"""
        self._errors.append(FastqError(name, error_type, error_message))

    def get_errors(self) -> List[FastqError]:
        """获取所有错误"""
        return self._errors


def extract_lib_id(lib_path: Path) -> str:
    name = lib_path.name
    lib_id = name.split("-")[-1]
    if lib_id.isdigit():
        lib_id = "-".join(name.split("-")[-2:])
    return lib_id


class FastqProcessor:
    """FASTQ文件处理器"""

    def __init__(
        self, base_dir: Path, error_recorder: FastqErrorRecorder
    ):  # 使用Optional类型提示
        self.base_dir = Path(base_dir)
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

        lib_id = extract_lib_id(sample_path)
        lib_info = []

        if not fastqs:
            logger.warning(f"在 {sample_path} 中未找到FASTQ文件")
            self.error_recorder.record_error(
                lib_id,
                FastqErrorType.INCOMPLETE.value,  # 使用value属性
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
                    error_type=FastqErrorType.FORMAT.value,  # 使用value属性
                    error_message=f"无法识别的FASTQ文件: {fq.name}",
                )

        r1_count = sum(1 for fq in fastqs if fq.name.endswith("_R1.fastq.gz"))
        r2_count = sum(1 for fq in fastqs if fq.name.endswith("_R2.fastq.gz"))

        if r1_count != r2_count:
            self.error_recorder.record_error(
                name=lib_id,
                error_type=FastqErrorType.FORMAT.value,
                error_message=f"R1和R2数量不一致: {r1_count} != {r2_count}",
            )

        return lib_info

    def _libid_not_duplicated(self, fq_line: Path) -> None:
        try:
            sample_dirs = list(fq_line.glob("Sample*"))
        except Exception as e:
            raise Exception(f"遍历目录失败: {e}")

        libids = [extract_lib_id(each_dir) for each_dir in sample_dirs]
        libid_counts = Counter(libids)
        duplicated = {libid for libid, count in libid_counts.items() if count > 1}

        if duplicated:
            # raise Exception(f"发现重复的 libid: {duplicated}，请检查目录命名！")
            for lib_id in duplicated:
                self.error_recorder.record_error(
                    name=lib_id,
                    error_type=FastqErrorType.DUPLICATED.value,
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

    def read_or_build_config(
        self, fq_line_dir: Path, force_rebuild: bool = False
    ) -> pd.DataFrame:
        """读取或构建配置文件"""
        config_file = fq_line_dir / "libid_fastq_config.tsv"

        if not force_rebuild and config_file.exists():
            try:
                logger.info(f"使用已有配置文件: {config_file}")
                df = pd.read_table(config_file)
                # 验证必要的列是否存在
                required_cols = ["libid", "read_type", "path"]
                if not all(col in df.columns for col in required_cols):
                    logger.warning(f"配置文件缺少必要列，重新构建: {config_file}")
                    raise ValueError("配置文件格式不正确")

                check_lib_map(self.error_recorder, df, fq_line_dir)
                return df
            except Exception as e:
                logger.warning(f"读取配置文件失败，重新构建: {e}")

        libid_map = self.build_libid_fastq_map(fq_line_dir)
        check_lib_map(self.error_recorder, libid_map, fq_line_dir)

        if not libid_map.empty:
            libid_map["dir_name"] = fq_line_dir.name
            try:
                libid_map.to_csv(config_file, sep="\t", index=False)
                logger.success(f"配置文件已保存: {config_file}")
            except Exception as e:
                logger.error(f"保存配置文件失败: {e}")

        return libid_map

    def load_config(
        self, fq_lines: np.ndarray, force_rebuild: bool = False
    ) -> pd.DataFrame:
        """加载所有相关的配置"""
        libid_map_list = []
        target_dirs = []

        # 查找所有匹配的目录
        for date_dir in self.base_dir.glob("20*"):
            if not date_dir.is_dir():
                continue
            for tcwl_dir in date_dir.glob("tcwl-*"):
                if tcwl_dir.name in fq_lines:
                    target_dirs.append(tcwl_dir)

        if not target_dirs:
            logger.error("未找到匹配的目录")
            return pd.DataFrame()

        logger.info(f"找到 {len(target_dirs)} 个匹配目录")

        for each_path in tqdm(target_dirs, desc="加载配置"):
            # self._libid_not_duplicated(each_path)
            try:
                logger.info(f"获取libid-fastq配置：{each_path.name}")
                libid_map = self.read_or_build_config(each_path, force_rebuild)
                if not libid_map.empty:
                    libid_map_list.append(libid_map)
            except Exception as e:
                logger.error(f"处理 {each_path} 时出错: {e}")
                continue

        if not libid_map_list:
            logger.error("未获取到任何有效配置")
            return pd.DataFrame()

        all_libid_map = pd.concat(libid_map_list, ignore_index=True)
        logger.success(f"成功加载 {len(all_libid_map)} 条配置记录")
        return all_libid_map


class ScriptRunner:
    """脚本执行器"""

    @staticmethod
    def merge_or_link_command(fq_list: List[str], output_name: str) -> str:
        """生成合并或链接命令"""
        if len(fq_list) == 1:
            return f"cp {fq_list[0]} {output_name}"
        return f"cat {' '.join(fq_list)} > {output_name}"

    @staticmethod
    def run_script(script_path: Path) -> Tuple[bool, str]:
        """运行单个脚本"""
        try:
            subprocess.run(
                ["bash", str(script_path)], check=True, capture_output=True, text=True
            )
            return True, f"成功: {script_path.name}"
        except subprocess.CalledProcessError as e:
            error_msg = f"失败: {script_path.name} - {e.stderr}"
            logger.error(error_msg)
            return False, error_msg

    @staticmethod
    def run_scripts_in_parallel(
        scripts_dir: Path, max_workers: int = 8
    ) -> Dict[str, int]:
        """并行运行脚本"""
        scripts = list(scripts_dir.glob("mergeFastq-*.sh"))
        if not scripts:
            logger.warning(f"在 {scripts_dir} 中未找到脚本文件")
            return {"success": 0, "failed": 0}

        logger.info(f"准备并行执行 {len(scripts)} 个脚本")

        results = {"success": 0, "failed": 0}

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(ScriptRunner.run_script, script): script
                for script in scripts
            }

            for future in tqdm(
                as_completed(futures), total=len(futures), desc="执行脚本"
            ):
                success, message = future.result()
                if success:
                    results["success"] += 1
                else:
                    results["failed"] += 1

        logger.info(
            f"脚本执行完成: 成功 {results['success']}, 失败 {results['failed']}"
        )
        return results


def write_nextflow_input(
    fq_df: pd.DataFrame,
    output_dir: Path,
    error_recorder: FastqErrorRecorder,
    warning_recorder: FastqErrorRecorder,
    threads: int = 8,
    run_script: bool = True,
) -> Optional[Dict[str, int]]:
    """写入Nextflow输入文件"""
    if fq_df.empty:
        logger.error("没有数据可写入")
        return None

    scripts_dir = output_dir / "scripts"
    scripts_dir.mkdir(exist_ok=True, parents=True)

    script_count = 0
    miss_df = fq_df[fq_df["path"].isna()]
    if not miss_df.empty:
        for row in miss_df.itertuples():
            error_recorder.record_error(
                name=str(row.libid),
                error_type="INCOMPLETE",  # 直接使用字符串
                error_message=f"{row.libid}: {row.libid}-{row.sample_id}-{row.dir_name} 没有找到数据",
            )

    for (sample_id, read_type), sample_df in fq_df.groupby(["sample_id", "read_type"]):
        miss_df = sample_df[sample_df["path"].isna()]
        if not miss_df.empty:
            logger.warning(f"包含缺失路径的样品: {sample_id}-{read_type}")
            error_recorder.record_error(
                name=sample_id,
                error_type=FastqErrorType.INCOMPLETE.value,
                error_message=f"包含缺失路径的样品: {sample_id}-{read_type}",
            )
            continue

        out_fq = output_dir.absolute() / f"{sample_id}.{read_type}.fq.gz"
        fq_list = sorted(sample_df["path"].tolist())

        cmd_file = scripts_dir / f"mergeFastq-{sample_id}-{read_type}.sh"
        cmd = ScriptRunner.merge_or_link_command(fq_list, str(out_fq))

        try:
            with open(cmd_file, "w") as f:
                f.write(f"#!/bin/bash\n")
                f.write(f"set -euo pipefail\n")
                f.write(f"{cmd}\n")
            cmd_file.chmod(0o755)
            script_count += 1
        except Exception as e:
            logger.error(f"写入脚本文件失败 {cmd_file}: {e}")

    logger.info(f"生成了 {script_count} 个脚本文件")

    errors = error_recorder.get_errors()
    if errors:
        for each_error in errors:
            logger.error(f"{each_error.error_type} - {each_error.error_message}")
    warnings = warning_recorder.get_errors()
    if warnings:
        for each_warning in warnings:
            logger.warning(f"{each_warning.error_type} - {each_warning.error_message}")

    should_run_script = run_script and script_count > 0 and len(errors) == 0

    if should_run_script:
        return ScriptRunner.run_scripts_in_parallel(scripts_dir, max_workers=threads)
    else:
        logger.warning("发现错误，跳过脚本执行")

    return None


def validate_sample_info(sample_df: pd.DataFrame) -> pd.DataFrame:
    """验证样品信息"""
    required_cols = ["libid", "sample_id", "dir_name"]
    missing_cols = [col for col in required_cols if col not in sample_df.columns]

    if missing_cols:
        raise ValueError(f"样品信息文件缺少必要列: {missing_cols}")

    # 去除重复行
    original_count = len(sample_df)
    sample_df = sample_df.drop_duplicates()
    if len(sample_df) < original_count:
        logger.warning(f"去除了 {original_count - len(sample_df)} 行重复数据")

    # 检查空值
    null_counts = sample_df[required_cols].isnull().sum()
    if null_counts.any():
        logger.warning(f"发现空值: {null_counts[null_counts > 0].to_dict()}")

    return sample_df


def log_statistics(df: pd.DataFrame):
    """记录统计信息"""
    if df.empty:
        logger.error("数据框为空")
        return

    total_samples = df["sample_id"].nunique()
    total_libs = df["libid"].nunique()

    miss_df = df[df["path"].isna()]
    miss_samples = miss_df["sample_id"].nunique()
    miss_libs = miss_df["libid"].nunique()

    logger.info(f"总计: {total_samples} 个样品, {total_libs} 个文库")

    if miss_samples > 0:
        logger.error(f"缺失数据: {miss_samples} 个样品, {miss_libs} 个文库")
        missing_sample_ids = miss_df["sample_id"].unique()[:10]  # 只显示前10个
        logger.error(f"缺失样品示例: {missing_sample_ids}")
    else:
        logger.success("所有样品数据已找到")


def check_sample_map(
    error_recorder: FastqErrorRecorder,
    warning_recorder: FastqErrorRecorder,
    df: pd.DataFrame,
    low_data_threshold=0.01,
):
    logger.info("正在检查样本映射关系...")
    duplicated_lines = df[df.duplicated(subset=["libid", "sample_id", "dir_name"])]
    if not duplicated_lines.empty:
        logger.warning("样本映射关系有重复项:")
        for row in duplicated_lines.itertuples():
            warning_recorder.record_error(
                name=str(row.libid),
                error_type="DUPLICATED",  # 直接使用字符串
                error_message=f"样本映射关系有重复项: {row.libid}-{row.sample_id}-{row.dir_name}",
            )
    low_data_df = df[df["data_size"] < low_data_threshold]
    if not low_data_df.empty:
        logger.warning(f"{len(low_data_df)}个样本数据量小于{low_data_threshold}G")
        for row in low_data_df.itertuples():
            error_recorder.record_error(
                name=str(row.libid),
                error_type="INCOMPLETE",  # 直接使用字符串
                error_message=f"样本数据量小于{low_data_threshold}G: {row.libid}-{row.sample_id}-{row.dir_name}",
            )


def check_lib_map(
    error_reccoder: FastqErrorRecorder, df: pd.DataFrame, fq_line_dir: Path
):
    r1_r2_count_df = df[["libid", "read_type"]].value_counts().unstack(1).reset_index()
    r1_r2_ne_df = r1_r2_count_df[r1_r2_count_df["R1"] != r1_r2_count_df["R2"]]
    if not r1_r2_ne_df.empty:
        for row in r1_r2_ne_df.itertuples():
            error_reccoder.record_error(
                name=str(row.libid),
                error_type="INCOMPLETE",  # 直接使用字符串
                error_message=f"文库映射关系有误: {fq_line_dir} - {row.libid} - R1 R2 not equal",
            )


@app.command()
def run(
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

    try:
        # 验证输入文件
        if not sample_info.exists():
            logger.error(f"样品信息文件不存在: {sample_info}")
            raise typer.Exit(1)

        # 读取样品信息
        logger.info(f"读取样品信息: {sample_info}")
        try:
            sample_df = pd.read_table(
                sample_info,
                header=None,
                names=["sample_id", "data_size", "dir_name", "libid"],
                usecols=[1, 2, 3, 5],
            )
        except Exception as e:
            logger.error(f"读取样品信息文件失败: {e}")
            raise typer.Exit(1)

        # 验证样品信息
        sample_df = validate_sample_info(sample_df)
        sample_libs = sample_df["dir_name"].unique()
        logger.info(f"需要处理 {len(sample_libs)} 个数据目录")

        # 初始化错误收集器
        error_collector = FastqErrorRecorder()
        warning_collector = FastqErrorRecorder()
        processor = FastqProcessor(base_dir, error_recorder=error_collector)

        check_sample_map(
            error_collector, warning_collector, sample_df, empty_data_threshold
        )  # 使用模块级函数
        sample_df = sample_df.drop_duplicates(subset=["sample_id", "dir_name", "libid"])

        # 加载配置
        logger.info("加载FASTQ文件配置")
        libid_map = processor.load_config(sample_libs, force_rebuild=force_rebuild)

        if libid_map.empty:
            logger.error("未找到任何FASTQ文件配置")
            raise typer.Exit(1)

        # 合并数据
        logger.info("合并样品信息和FASTQ配置")
        merged_df = sample_df.merge(libid_map, how="left")

        # 记录统计信息
        log_statistics(merged_df)

        # 保存检查结果
        try:
            merged_df.to_csv(check_file, sep="\t", index=False)
            logger.success(f"检查结果已保存: {check_file}")
        except Exception as e:
            logger.error(f"保存检查文件失败: {e}")

        # 生成输出文件
        if output_dir is not None:
            output_dir.mkdir(exist_ok=True, parents=True)
            logger.info(f"生成Nextflow输入文件到: {output_dir}")

            results = write_nextflow_input(
                merged_df, output_dir, error_collector, threads=threads
            )

            if results:
                logger.info(f"脚本执行结果: {results}")

            logger.success(f"处理完成！输出目录: {output_dir}")

    except KeyboardInterrupt:
        logger.info("用户中断操作")
        raise typer.Exit(130)
    except Exception as e:
        logger.error(f"处理过程中出现错误: {e}")
        raise typer.Exit(1)


@app.command()
def validate(
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

    try:
        # 验证输入文件
        if not sample_info.exists():
            logger.error(f"样品信息文件不存在: {sample_info}")
            raise typer.Exit(1)

        # 读取样品信息
        logger.info(f"读取样品信息: {sample_info}")
        try:
            sample_df = pd.read_table(
                sample_info,
                header=None,
                names=["sample_id", "data_size", "dir_name", "libid"],
                usecols=[1, 2, 3, 5],
            )
        except Exception as e:
            logger.error(f"读取样品信息文件失败: {e}")
            raise typer.Exit(1)

        # 验证样品信息
        sample_df = validate_sample_info(sample_df)
        sample_libs = sample_df["dir_name"].unique()
        logger.info(f"需要处理 {len(sample_libs)} 个数据目录")

        # 初始化错误收集器
        error_collector = FastqErrorRecorder()
        warning_collector = FastqErrorRecorder()
        processor = FastqProcessor(base_dir, error_recorder=error_collector)

        check_sample_map(
            error_collector, warning_collector, sample_df, empty_data_threshold
        )  # 使用模块级函数

        sample_df = sample_df.drop_duplicates(subset=["sample_id", "dir_name", "libid"])

        # 加载配置
        logger.info("加载FASTQ文件配置")
        libid_map = processor.load_config(sample_libs, force_rebuild=force_rebuild)

        if libid_map.empty:
            logger.error("未找到任何FASTQ文件配置")
            raise typer.Exit(1)

        # 合并数据
        logger.info("合并样品信息和FASTQ配置")
        merged_df = sample_df.merge(libid_map, how="left")

        # 记录统计信息
        log_statistics(merged_df)

        # 保存检查结果
        try:
            merged_df.to_csv(check_file, sep="\t", index=False)
            logger.success(f"检查结果已保存: {check_file}")
        except Exception as e:
            logger.error(f"保存检查文件失败: {e}")

        # 生成输出文件
        if output_dir is not None:
            output_dir.mkdir(exist_ok=True, parents=True)
            logger.info(f"生成Nextflow输入文件到: {output_dir}")

            write_nextflow_input(
                merged_df,
                output_dir,
                error_collector,
                threads=threads,
                run_script=False,
            )

            errors = error_collector.get_errors()
            if errors:
                for each_error in errors:
                    logger.error(
                        f"{each_error.error_type} - {each_error.error_message}"
                    )
            warnings = warning_collector.get_errors()
            if warnings:
                for each_warning in warnings:
                    logger.warning(
                        f"{each_warning.error_type} - {each_warning.error_message}"
                    )
            else:
                logger.success(f"检查完成：没有发现问题！")

    except KeyboardInterrupt:
        logger.info("用户中断操作")
        raise typer.Exit(130)
    except Exception as e:
        logger.error(f"处理过程中出现错误: {e}")
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
