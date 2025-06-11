#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import typer
from loguru import logger
from tqdm import tqdm

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


class FastqProcessor:
    """FASTQ文件处理器"""

    def __init__(self, base_dir: Path):
        self.base_dir = Path(base_dir)
        if not self.base_dir.exists():
            raise FileNotFoundError(f"基础目录不存在: {base_dir}")

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

        if not fastqs:
            logger.warning(f"在 {sample_path} 中未找到FASTQ文件")
            return []

        lib_id = filename.split("-")[-1]
        lib_info = []

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

        return lib_info

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
        sample_dirs = list(fastq_path.glob("Sample*"))

        if not sample_dirs:
            logger.warning(f"在 {fastq_path} 中未找到Sample*目录")
            return pd.DataFrame()

        for path in tqdm(sample_dirs, desc=f"处理 {fastq_path.name}"):
            try:
                info = self.parse_fastq_filename(path)
                libid_map.extend(info)
            except Exception as e:
                logger.error(f"处理 {path} 时出错: {e}")
                continue

        return pd.DataFrame(libid_map)

    def read_or_build_config(self, fq_line_dir: Path) -> pd.DataFrame:
        """读取或构建配置文件"""
        config_file = fq_line_dir / "libid_fastq_config.tsv"

        if config_file.exists():
            try:
                logger.info(f"使用已有配置文件: {config_file}")
                df = pd.read_table(config_file)
                # 验证必要的列是否存在
                required_cols = ["libid", "read_type", "path"]
                if not all(col in df.columns for col in required_cols):
                    logger.warning(f"配置文件缺少必要列，重新构建: {config_file}")
                    raise ValueError("配置文件格式不正确")
                return df
            except Exception as e:
                logger.warning(f"读取配置文件失败，重新构建: {e}")

        libid_map = self.build_libid_fastq_map(fq_line_dir)
        if not libid_map.empty:
            libid_map["dir_name"] = fq_line_dir.name
            try:
                libid_map.to_csv(config_file, sep="\t", index=False)
                logger.success(f"配置文件已保存: {config_file}")
            except Exception as e:
                logger.error(f"保存配置文件失败: {e}")

        return libid_map

    def load_config(self, fq_lines: np.ndarray) -> pd.DataFrame:
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
            try:
                logger.info(f"获取libid-fastq配置：{each_path.name}")
                libid_map = self.read_or_build_config(each_path)
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
            return f"ln -sf {fq_list[0]} {output_name}"  # 使用软链接节省空间
        return f"cat {' '.join(fq_list)} > {output_name}"

    @staticmethod
    def run_script(script_path: Path) -> Tuple[bool, str]:
        """运行单个脚本"""
        try:
            result = subprocess.run(
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
    fq_df: pd.DataFrame, output_dir: Path, run: bool = False, threads: int = 8
) -> Optional[Dict[str, int]]:
    """写入Nextflow输入文件"""
    if fq_df.empty:
        logger.error("没有数据可写入")
        return None

    scripts_dir = output_dir / "scripts"
    scripts_dir.mkdir(exist_ok=True, parents=True)

    script_count = 0

    for (sample_id, read_type), sample_df in fq_df.groupby(["sample_id", "read_type"]):
        if sample_df["path"].isna().any():
            logger.warning(f"跳过缺失路径的样品: {sample_id}-{read_type}")
            continue

        out_fq = output_dir.absolute() / f"{sample_id}.{read_type}.fq.gz"
        fq_list = sorted(sample_df["path"].dropna().tolist())

        if not fq_list:
            logger.warning(f"样品 {sample_id}-{read_type} 没有有效的FASTQ文件")
            continue

        cmd_file = scripts_dir / f"mergeFastq-{sample_id}-{read_type}.sh"
        cmd = ScriptRunner.merge_or_link_command(fq_list, str(out_fq))

        try:
            with open(cmd_file, "w") as f:
                f.write(f"#!/bin/bash\n")
                f.write(f"set -euo pipefail\n")  # 严格错误处理
                f.write(f"{cmd}\n")
            cmd_file.chmod(0o755)  # 使脚本可执行
            script_count += 1
        except Exception as e:
            logger.error(f"写入脚本文件失败 {cmd_file}: {e}")

    logger.info(f"生成了 {script_count} 个脚本文件")

    if run and script_count > 0:
        return ScriptRunner.run_scripts_in_parallel(scripts_dir, max_workers=threads)

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


@app.command()
def main(
    sample_info: Path = typer.Argument(
        ..., help="样品信息TSV文件，必须包含libid、sample_id、dir_name列"
    ),
    base_dir: Path = typer.Option(DEFAULT_BASE_DIR, help="包含所有FASTQ数据的基础目录"),
    output_dir: Optional[Path] = typer.Option(
        Path("raw_data"), help="FASTQ文件输出目录"
    ),
    check_file: Path = typer.Option("check_file.tsv", help="检查结果输出文件"),
    threads: int = typer.Option(8, min=1, max=32, help="并行处理线程数"),
    run: bool = typer.Option(False, help="是否立即执行合并脚本"),
    force_rebuild: bool = typer.Option(False, help="强制重建配置文件"),
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
                names=["libid", "sample_id", "data_size", "dir_name"],
                usecols=[0, 1, 2, 3],
            )
        except Exception as e:
            logger.error(f"读取样品信息文件失败: {e}")
            raise typer.Exit(1)

        # 验证样品信息
        sample_df = validate_sample_info(sample_df)
        sample_libs = sample_df["dir_name"].unique()
        logger.info(f"需要处理 {len(sample_libs)} 个数据目录")

        # 初始化处理器
        processor = FastqProcessor(base_dir)

        # 如果需要强制重建，删除现有配置文件
        if force_rebuild:
            logger.info("强制重建模式：删除现有配置文件")
            for date_dir in base_dir.glob("20*"):
                for tcwl_dir in date_dir.glob("tcwl-*"):
                    config_file = tcwl_dir / "libid_fastq_config.tsv"
                    if config_file.exists():
                        config_file.unlink()

        # 加载配置
        logger.info("加载FASTQ文件配置")
        libid_map = processor.load_config(sample_libs)

        if libid_map.empty:
            logger.error("未找到任何FASTQ文件配置")
            raise typer.Exit(1)

        # 合并数据
        logger.info("合并样品信息和FASTQ配置")
        merged_df = sample_df.merge(libid_map, on="libid", how="left")

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
                merged_df, output_dir, run=run, threads=threads
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
    sample_info: Path = typer.Argument(..., help="样品信息TSV文件"),
    base_dir: Path = typer.Option(DEFAULT_BASE_DIR, help="基础数据目录"),
):
    """验证样品信息和数据完整性"""

    if not sample_info.exists():
        logger.error(f"文件不存在: {sample_info}")
        raise typer.Exit(1)

    try:
        sample_df = pd.read_table(
            sample_info,
            header=None,
            names=["libid", "sample_id", "data_size", "dir_name"],
            usecols=[0, 1, 2, 3],
        )
        sample_df = validate_sample_info(sample_df)

        processor = FastqProcessor(base_dir)
        libid_map = processor.load_config(sample_df["dir_name"].unique())

        merged_df = sample_df.merge(libid_map, on="libid", how="left")
        log_statistics(merged_df)

        logger.success("验证完成")

    except Exception as e:
        logger.error(f"验证失败: {e}")
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
