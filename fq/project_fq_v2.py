#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""FASTQ文件处理和合并工具 v2.0"""

import subprocess
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Annotated

import numpy as np
import pandas as pd
import typer
from loguru import logger
from tqdm import tqdm

__version__ = "2.0"

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


class IncludeExcludePathOverlapError(Exception):
    """包含和排除路径重叠错误"""


class FastqErrorType(StrEnum):
    DUPLICATED = "DUPLICATED"
    INCOMPLETE = "INCOMPLETE"
    FORMAT = "FORMAT"
    DATA_SIZE = "DATA_SIZE"
    LOW_DATA = "LOW_DATA"


class DataMode(StrEnum):
    cp = "cp"
    link = "link"


@dataclass
class FastqError:
    name: str
    error_type: FastqErrorType
    error_message: str


@dataclass
class FastqErrorRecorder:
    _errors: list[FastqError] = field(default_factory=list)

    def record(self, name: str, error_type: FastqErrorType, message: str):
        self._errors.append(FastqError(name, error_type, message))

    @property
    def errors(self) -> list[FastqError]:
        return self._errors

    def __len__(self) -> int:
        return len(self._errors)

    def __bool__(self) -> bool:
        return len(self._errors) > 0


def extract_lib_id(lib_path: Path) -> str:
    name = lib_path.name
    lib_id = name.split("-")[-1]
    if lib_id.isdigit() or (len(lib_id) == 1 and lib_id.islower()):
        lib_id = "-".join(name.split("-")[-2:])
    return lib_id


def load_path_set(file_path: Path | None) -> set[Path]:
    """将文件中的路径读取为Path集合"""
    if file_path is None:
        return set()
    if not file_path.exists():
        raise FileNotFoundError(f"-i/-e 文件列表文件不存在: {file_path}")

    path_set = set()
    for line in file_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line:
            line_path = Path(line).resolve()
            if not line_path.exists():
                raise FileNotFoundError(f"-i/-e 列表中的路径不存在: {line_path}")
            path_set.add(line_path)
    return path_set


def is_path_included(path: Path, path_set: set[Path]) -> bool:
    """检查路径是否在集合中或是其子路径"""
    resolved = path.resolve()
    if resolved in path_set:
        return True
    return any(resolved.is_relative_to(p) for p in path_set)


class FastqProcessor:
    """FASTQ文件处理器"""

    def __init__(
        self,
        base_dir: Path,
        error_recorder: FastqErrorRecorder,
        exclude_paths: Path | None = None,
        include_paths: Path | None = None,
    ):
        self.base_dir = Path(base_dir)
        self.error_recorder = error_recorder
        self.line_tracker: list[dict] = []
        self.duplicated_data_df = pd.DataFrame()
        self.include_path_set = load_path_set(include_paths)
        self.exclude_path_set = load_path_set(exclude_paths)
        self._validate_paths()

    def _validate_paths(self) -> None:
        overlap = self.include_path_set & self.exclude_path_set
        if overlap:
            raise IncludeExcludePathOverlapError(
                f"包含和排除路径重叠: {', '.join(str(p) for p in overlap)}"
            )

    def _determine_read_type(self, filename: str) -> str | None:
        for pattern, read_type in READ_TYPE_PATTERNS.items():
            if filename.endswith(pattern):
                return read_type
        return None

    def parse_fastq_filename(self, sample_path: Path) -> list[dict]:
        """解析FASTQ文件名"""
        if not sample_path.exists():
            logger.warning(f"样品路径不存在: {sample_path}")
            return []

        fastqs = []
        for pattern in FASTQ_EXTENSIONS:
            fastqs.extend(sample_path.glob(pattern))

        lib_id = extract_lib_id(sample_path)

        if not fastqs:
            logger.warning(f"在 {sample_path} 中未找到FASTQ文件")
            self.error_recorder.record(
                lib_id, FastqErrorType.INCOMPLETE, f"未找到FASTQ文件 {sample_path}"
            )
            return []

        lib_info = []
        for fq in fastqs:
            read_type = self._determine_read_type(fq.name)
            if read_type:
                lib_info.append({
                    "libid": lib_id,
                    "read_type": read_type,
                    "path": str(fq.absolute()),
                })
            else:
                logger.warning(f"无法识别的FASTQ文件: {fq.name}")
                self.error_recorder.record(
                    lib_id, FastqErrorType.FORMAT, f"无法识别的FASTQ文件: {fq.name}"
                )

        # 检查R1/R2配对
        r1_count = sum(1 for fq in fastqs if "_R1.fastq.gz" in fq.name)
        r2_count = sum(1 for fq in fastqs if "_R2.fastq.gz" in fq.name)
        if r1_count != r2_count:
            self.error_recorder.record(
                lib_id, FastqErrorType.FORMAT, f"R1和R2数量不一致: {r1_count} != {r2_count}"
            )

        return lib_info

    def build_libid_fastq_map(self, fastq_path: Path) -> pd.DataFrame:
        """构建library ID到FASTQ文件的映射"""
        if not fastq_path.exists():
            logger.error(f"FASTQ路径不存在: {fastq_path}")
            return pd.DataFrame()

        sample_dirs = list(fastq_path.glob("Sample*"))
        if not sample_dirs:
            raise ValueError(f"在 {fastq_path} 中未找到Sample*目录")

        libid_map = []
        for path in tqdm(sample_dirs, desc=f"处理 {fastq_path.name}"):
            info = self.parse_fastq_filename(path)
            libid_map.extend(info)

        return pd.DataFrame(libid_map)

    def track_line_lib_directory(self, fastq_path: Path) -> None:
        sample_dirs = list(fastq_path.glob("Sample*"))
        for each_path in sample_dirs:
            self.line_tracker.append({
                "line": fastq_path.name,
                "lib_dir": each_path.name,
                "line_path": str(fastq_path.absolute()),
            })

    def check_duplicated_data(self) -> None:
        """检查重复数据"""
        line_track_df = pd.DataFrame(self.line_tracker)
        dup_df = line_track_df[
            line_track_df.duplicated(subset=["line", "lib_dir"], keep=False)
        ]

        if dup_df.empty:
            return

        for line, df_j in dup_df.groupby("line"):
            if len(df_j) > 1:
                dup_lib_dirs = df_j["lib_dir"].unique().tolist()
                dir_names = ",".join(dup_lib_dirs[:3])
                if len(dup_lib_dirs) > 3:
                    dir_names = f"{dir_names} ...，共{len(dup_lib_dirs)}个"
                dup_paths = " | ".join(df_j["line_path"].unique().tolist())
                self.error_recorder.record(
                    str(line),
                    FastqErrorType.DUPLICATED,
                    f"<cyan>文库目录</cyan> {dup_paths} <cyan>包含重复的数据:</cyan> <w>{dir_names}</w>",
                )

        self.duplicated_data_df = (
            dup_df.groupby(["line", "lib_dir"])["line_path"]
            .unique()
            .map(" | ".join)
            .reset_index()
        )

    def read_or_build_config(self, fq_line_dir: Path, force_rebuild: bool = False) -> pd.DataFrame:
        """读取或构建配置文件"""
        config_file = fq_line_dir / "libid_fastq_config.tsv"
        self.track_line_lib_directory(fq_line_dir)

        if not force_rebuild and config_file.exists():
            try:
                logger.info(f"使用已有配置文件: {config_file}")
                df = pd.read_table(config_file)
                required_cols = ["libid", "read_type", "path"]
                if all(col in df.columns for col in required_cols):
                    check_lib_map(self.error_recorder, df, fq_line_dir)
                    return df
                logger.warning(f"配置文件缺少必要列，重新构建: {config_file}")
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

    def load_config(self, fq_lines: np.ndarray, force_rebuild: bool = False) -> pd.DataFrame:
        """加载所有相关的配置"""
        target_dirs = []

        # 查找所有匹配的目录
        for date_dir in self.base_dir.glob("20*"):
            if not date_dir.is_dir():
                continue
            for tcwl_dir in date_dir.glob("*"):
                if tcwl_dir.name in fq_lines and not is_path_included(tcwl_dir, self.exclude_path_set):
                    target_dirs.append(tcwl_dir.resolve())

        # 添加include路径
        for include_path in self.include_path_set:
            for child in include_path.glob("*"):
                if child.is_dir() and child.resolve() not in target_dirs:
                    target_dirs.append(child.resolve())

        if not target_dirs:
            logger.error("未找到匹配的目录")
            return pd.DataFrame()

        logger.info(f"找到 {len(target_dirs)} 个匹配目录")

        libid_map_list = []
        for each_path in tqdm(target_dirs, desc="加载配置"):
            try:
                logger.info(f"获取libid-fastq配置：{each_path.name}")
                libid_map = self.read_or_build_config(each_path, force_rebuild)
                if not libid_map.empty:
                    libid_map_list.append(libid_map)
            except Exception as e:
                logger.error(f"处理 {each_path} 时出错: {e}")

        self.check_duplicated_data()

        if not libid_map_list:
            logger.error("未获取到任何有效配置")
            return pd.DataFrame()

        all_libid_map = pd.concat(libid_map_list, ignore_index=True)
        logger.success(f"成功加载 {len(all_libid_map)} 条配置记录")
        return all_libid_map


class ScriptRunner:
    """脚本执行器"""

    @staticmethod
    def merge_or_link_command(fq_list: list[str], output_name: str, mode: DataMode) -> str:
        if len(fq_list) == 1:
            cmd = "ln -s" if mode == DataMode.link else "cp"
            return f"{cmd} {fq_list[0]} {output_name}"
        return f"cat {' '.join(fq_list)} > {output_name}"

    @staticmethod
    def run_script(script_path: Path) -> tuple[bool, str]:
        try:
            subprocess.run(["bash", str(script_path)], check=True, capture_output=True, text=True)
            return True, f"成功: {script_path.name}"
        except subprocess.CalledProcessError as e:
            error_msg = f"失败: {script_path.name} - {e.stderr}"
            logger.error(error_msg)
            return False, error_msg

    @staticmethod
    def run_scripts_in_parallel(scripts_dir: Path, max_workers: int = 8) -> dict[str, int]:
        scripts = list(scripts_dir.glob("mergeFastq-*.sh"))
        if not scripts:
            logger.warning(f"在 {scripts_dir} 中未找到脚本文件")
            return {"success": 0, "failed": 0}

        logger.info(f"准备并行执行 {len(scripts)} 个脚本")
        results = {"success": 0, "failed": 0}

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(ScriptRunner.run_script, s): s for s in scripts}
            for future in tqdm(as_completed(futures), total=len(futures), desc="执行脚本"):
                success, _ = future.result()
                results["success" if success else "failed"] += 1

        logger.info(f"脚本执行结果: 成功 {results['success']}, 失败 {results['failed']}")
        return results


def check_lib_map(error_recorder: FastqErrorRecorder, df: pd.DataFrame, fq_line_dir: Path):
    """检查文库映射的R1/R2配对"""
    if df.empty:
        return
    r1_r2_count_df = df[["libid", "read_type"]].value_counts().unstack(fill_value=0).reset_index()
    if "R1" not in r1_r2_count_df.columns or "R2" not in r1_r2_count_df.columns:
        return
    r1_r2_ne_df = r1_r2_count_df[r1_r2_count_df["R1"] != r1_r2_count_df["R2"]]
    for row in r1_r2_ne_df.itertuples():
        error_recorder.record(
            str(row.libid),
            FastqErrorType.INCOMPLETE,
            f"文库映射关系有误: {fq_line_dir} - {row.libid} - R1 R2 not equal",
        )


def check_sample_map(
    error_recorder: FastqErrorRecorder,
    warning_recorder: FastqErrorRecorder,
    df: pd.DataFrame,
    low_data_threshold: float = 0.01,
):
    """检查样本映射关系"""
    logger.info("正在检查样本映射关系...")

    # 检查重复项
    dup_mask = df.duplicated(subset=["libid", "sample_id", "dir_name"])
    for row in df[dup_mask].itertuples():
        warning_recorder.record(
            str(row.libid),
            FastqErrorType.DUPLICATED,
            f"样本映射关系有重复项: {row.libid} | {row.sample_id} | {row.dir_name}",
        )

    # 检查低数据量
    low_data_df = df[df["data_size"] < low_data_threshold]
    if not low_data_df.empty:
        logger.warning(f"{len(low_data_df)}个样本数据量小于{low_data_threshold}G")
        for row in low_data_df.itertuples():
            error_recorder.record(
                str(row.libid),
                FastqErrorType.INCOMPLETE,
                f"样本数据量小于{low_data_threshold}G: {row.libid} | {row.sample_id} | {row.dir_name}",
            )


def validate_sample_info(sample_df: pd.DataFrame) -> pd.DataFrame:
    """验证样品信息"""
    required_cols = ["libid", "sample_id", "dir_name"]
    missing_cols = [col for col in required_cols if col not in sample_df.columns]
    if missing_cols:
        raise ValueError(f"样品信息文件缺少必要列: {missing_cols}")

    original_count = len(sample_df)
    sample_df = sample_df.drop_duplicates()
    if len(sample_df) < original_count:
        logger.warning(f"去除了 {original_count - len(sample_df)} 行重复数据")

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

    logger.info(f"总计: {total_samples} 个样品, {total_libs} 个文库")

    if not miss_df.empty:
        logger.error(f"缺失数据: {miss_df['sample_id'].nunique()} 个样品, {miss_df['libid'].nunique()} 个文库")
        logger.error(f"缺失样品示例: {miss_df['sample_id'].unique()[:10]}")
    else:
        logger.success("所有样品数据已找到")


def write_nextflow_input(
    fq_df: pd.DataFrame,
    output_dir: Path,
    error_recorder: FastqErrorRecorder,
    warning_recorder: FastqErrorRecorder,
    threads: int = 8,
    run_script: bool = True,
    mode: DataMode = DataMode.link,
) -> dict[str, int] | None:
    """写入Nextflow输入文件"""
    if fq_df.empty:
        logger.error("没有数据可写入")
        return None

    scripts_dir = output_dir / "scripts"
    scripts_dir.mkdir(exist_ok=True, parents=True)

    # 记录缺失数据
    miss_df = fq_df[fq_df["path"].isna()]
    for row in miss_df.itertuples():
        error_recorder.record(
            str(row.libid),
            FastqErrorType.INCOMPLETE,
            f"{row.libid}: {row.libid}-{row.sample_id}-{row.dir_name} 没有找到数据",
        )

    script_count = 0
    for (sample_id, read_type), sample_df in fq_df.groupby(["sample_id", "read_type"]):
        if sample_df["path"].isna().any():
            logger.warning(f"包含缺失路径的样品: {sample_id}-{read_type}")
            error_recorder.record(
                str(sample_id),
                FastqErrorType.INCOMPLETE,
                f"包含缺失路径的样品: {sample_id}-{read_type}",
            )
            continue

        if run_script:
            out_fq = output_dir.absolute() / f"{sample_id}.{read_type}.fq.gz"
            fq_list = sorted(sample_df["path"].tolist())
            cmd = ScriptRunner.merge_or_link_command(fq_list, str(out_fq), mode)

            cmd_file = scripts_dir / f"mergeFastq-{sample_id}-{read_type}.sh"
            cmd_file.write_text(f"#!/bin/bash\nset -euo pipefail\n{cmd}\n", encoding="utf-8")
            cmd_file.chmod(0o755)
            script_count += 1

    if script_count:
        logger.info(f"生成了 {script_count} 个脚本文件")

    # 输出错误和警告
    for err in error_recorder.errors:
        logger.opt(colors=True).error(f"{err.error_type} - {err.error_message}")
    for warn in warning_recorder.errors:
        logger.warning(f"{warn.error_type} - {warn.error_message}")

    if run_script:
        if error_recorder:
            logger.warning("发现错误，跳过脚本执行")
            return None
        if script_count == 0:
            logger.warning("没有生成脚本，跳过脚本执行")
            return None
        return ScriptRunner.run_scripts_in_parallel(scripts_dir, max_workers=threads)

    if not error_recorder and not warning_recorder:
        logger.success("检查完成：没有发现问题！")
    return None


def _process_common(
    sample_info: Path,
    base_dir: Path,
    output_dir: Path,
    check_file: Path,
    threads: int,
    force_rebuild: bool,
    empty_data_threshold: float,
    exclude: Path | None,
    include: Path | None,
    mode: DataMode,
    run_script: bool,
) -> None:
    """run 和 validate 的公共处理逻辑"""
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

    sample_df = validate_sample_info(sample_df)
    sample_libs = sample_df["dir_name"].unique()
    logger.info(f"需要处理 {len(sample_libs)} 个数据目录")

    # 初始化错误收集器
    error_collector = FastqErrorRecorder()
    warning_collector = FastqErrorRecorder()
    processor = FastqProcessor(
        base_dir,
        error_recorder=error_collector,
        exclude_paths=exclude,
        include_paths=include,
    )

    check_sample_map(error_collector, warning_collector, sample_df, empty_data_threshold)
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
    log_statistics(merged_df)

    # 保存检查结果
    try:
        merged_df.to_csv(check_file, sep="\t", index=False)
        logger.success(f"检查结果已保存: {check_file}")
    except Exception as e:
        logger.error(f"保存检查文件失败: {e}")

    # 生成输出文件
    output_dir.mkdir(exist_ok=True, parents=True)
    logger.info(f"生成Nextflow输入文件到: {output_dir}")

    results = write_nextflow_input(
        merged_df,
        output_dir,
        error_collector,
        warning_collector,
        threads=threads,
        run_script=run_script,
        mode=mode,
    )

    if not processor.duplicated_data_df.empty:
        dup_file = output_dir / "duplicated_data.tsv"
        try:
            processor.duplicated_data_df.to_csv(dup_file, sep="\t", index=False)
            logger.success(f"重复数据详情已保存: {dup_file}")
        except Exception as e:
            logger.error(f"保存重复数据文件失败: {e}")

    if results:
        logger.info(f"脚本执行结果: {results}")

    logger.success(f"处理完成！输出目录: {output_dir}")


# 定义公共参数类型
SampleInfoArg = Annotated[Path, typer.Argument(help="样品信息TSV文件，必须包含libid、sample_id、dir_name列")]
BaseDirOpt = Annotated[Path, typer.Option(help="包含所有FASTQ数据的基础目录")]
OutputDirOpt = Annotated[Path, typer.Option(help="FASTQ文件输出目录")]
CheckFileOpt = Annotated[Path, typer.Option(help="检查结果输出文件")]
ThreadsOpt = Annotated[int, typer.Option(min=1, max=32, help="并行处理线程数")]
ForceRebuildOpt = Annotated[bool, typer.Option(help="强制重建配置文件")]
ThresholdOpt = Annotated[float, typer.Option(help="空数据阈值")]
ExcludeOpt = Annotated[Path | None, typer.Option("-e", "--exclude", help="需要排除的line路径文件")]
IncludeOpt = Annotated[Path | None, typer.Option("-i", "--include", help="需要包含的line路径文件")]


@app.command()
def run(
    sample_info: SampleInfoArg,
    base_dir: BaseDirOpt = DEFAULT_BASE_DIR,
    output_dir: OutputDirOpt = Path("."),
    check_file: CheckFileOpt = Path("check_file.tsv"),
    threads: ThreadsOpt = 8,
    force_rebuild: ForceRebuildOpt = False,
    empty_data_threshold: ThresholdOpt = 0.01,
    exclude: ExcludeOpt = None,
    include: IncludeOpt = None,
    mode: DataMode = DataMode.link,
):
    """处理FASTQ文件并执行合并脚本"""
    try:
        _process_common(
            sample_info, base_dir, output_dir, check_file, threads,
            force_rebuild, empty_data_threshold, exclude, include, mode,
            run_script=True,
        )
    except KeyboardInterrupt:
        logger.info("用户中断操作")
        raise typer.Exit(130)
    except typer.Exit:
        raise
    except Exception as e:
        logger.error(f"处理过程中出现错误: {e}")
        raise typer.Exit(1)


@app.command()
def validate(
    sample_info: SampleInfoArg,
    base_dir: BaseDirOpt = DEFAULT_BASE_DIR,
    output_dir: OutputDirOpt = Path("."),
    check_file: CheckFileOpt = Path("check_file.tsv"),
    threads: ThreadsOpt = 8,
    force_rebuild: ForceRebuildOpt = False,
    empty_data_threshold: ThresholdOpt = 0.01,
    exclude: ExcludeOpt = None,
    include: IncludeOpt = None,
):
    """仅验证FASTQ文件，不执行合并脚本"""
    try:
        _process_common(
            sample_info, base_dir, output_dir, check_file, threads,
            force_rebuild, empty_data_threshold, exclude, include, DataMode.link,
            run_script=False,
        )
    except KeyboardInterrupt:
        logger.info("用户中断操作")
        raise typer.Exit(130)
    except typer.Exit:
        raise
    except Exception as e:
        logger.error(f"处理过程中出现错误: {e}")
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
