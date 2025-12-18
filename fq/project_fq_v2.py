#!/usr/bin/env python3
"""
FASTQ文件处理和合并工具 v2.0

功能：根据下机数据目录名称和lib-id与样品名称的对应表格，
去下机数据目录中链接/合并下机数据，并将libid改为样品名称。

模块结构：
    - Config: 配置管理
    - Models: 数据模型 (FastqInfo, LineTrack, FastqError)
    - ErrorRecorder: 错误记录器
    - PathResolver: 路径解析器
    - FastqScanner: FASTQ文件扫描器
    - ConfigBuilder: 配置文件构建器
    - SampleValidator: 样本验证器
    - ScriptGenerator: 脚本生成器
    - ScriptRunner: 脚本执行器
    - Pipeline: 处理流水线
    - CLI: 命令行接口
"""

from __future__ import annotations

import subprocess
from collections.abc import Iterable, Iterator
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from enum import StrEnum, auto
from itertools import chain
from pathlib import Path
from typing import Annotated, NamedTuple, Self, cast

import pandas as pd
import typer
from loguru import logger
from tqdm import tqdm

__version__ = "2.0"


# ============================================================================
# 配置管理 (Config)
# ============================================================================


@dataclass(frozen=True)
class Config:
    """
    全局配置类 - 集中管理所有可配置参数

    遵循单一职责原则：只负责配置参数的存储和访问
    """

    # 路径配置
    base_dir: Path = Path("/public/home/zxchen/data_trans")
    log_file: Path = Path("project_fq.log")

    # FASTQ文件匹配模式
    fastq_extensions: tuple[str, ...] = ("*.fastq.gz", "*.fq.gz")

    # 读取类型模式映射
    read_type_patterns: dict[str, str] = field(
        default_factory=lambda: {
            "combined_R1.fastq.gz": "R1",
            "combined_R2.fastq.gz": "R2",
            "_R1.fastq.gz": "R1",
            "_R2.fastq.gz": "R2",
            "_1.fastq.gz": "R1",
            "_2.fastq.gz": "R2",
        }
    )

    # 必需的列名
    required_sample_cols: frozenset[str] = frozenset({"libid", "sample_id", "dir_name"})
    required_config_cols: frozenset[str] = frozenset({"libid", "read_type", "path"})

    # 样本信息文件列映射
    sample_info_columns: dict[str, int] = field(
        default_factory=lambda: {
            "sample_id": 1,
            "data_size": 2,
            "dir_name": 3,
            "libid": 5,
        }
    )

    # 默认参数
    default_threads: int = 8
    default_threshold: float = 0.01
    max_threads: int = 32

    # 文件名模式
    config_filename: str = "libid_fastq_config.tsv"
    script_prefix: str = "mergeFastq"
    date_dir_pattern: str = "20*"
    sample_dir_pattern: str = "Sample*"


# 全局默认配置
DEFAULT_CONFIG = Config()


# ============================================================================
# 枚举类型 (Enums)
# ============================================================================


class FastqErrorType(StrEnum):
    """FASTQ文件错误类型"""

    DUPLICATED = auto()
    INCOMPLETE = auto()
    FORMAT = auto()
    DATA_SIZE = auto()
    LOW_DATA = auto()


class DataMode(StrEnum):
    """数据处理模式"""

    CP = "cp"
    LINK = "link"


# ============================================================================
# 异常类 (Exceptions)
# ============================================================================


class PathOverlapError(Exception):
    """包含和排除路径重叠错误"""


class ConfigurationError(Exception):
    """配置错误"""


class ValidationError(Exception):
    """验证错误"""


# ============================================================================
# 数据模型 (Models)
# ============================================================================


class FastqInfo(NamedTuple):
    """FASTQ文件信息 - 不可变数据结构"""

    libid: str
    read_type: str
    path: str


class LineTrack(NamedTuple):
    """目录跟踪信息 - 不可变数据结构"""

    line: str
    lib_dir: str
    line_path: str


@dataclass(frozen=True, slots=True)
class FastqError:
    """FASTQ文件错误记录 - 不可变数据结构"""

    name: str
    error_type: FastqErrorType
    message: str

    def log(self, level: str = "error") -> None:
        """输出日志"""
        log_func = getattr(logger.opt(colors=True), level)
        log_func(f"{self.error_type} - {self.message}")


# ============================================================================
# 错误记录器 (ErrorRecorder)
# ============================================================================


@dataclass
class ErrorRecorder:
    """
    错误记录器 - 收集和管理错误信息

    遵循单一职责原则：只负责错误的记录和输出
    支持迭代协议和布尔判断
    """

    _items: list[FastqError] = field(default_factory=list)

    def record(self, name: str, error_type: FastqErrorType, message: str) -> Self:
        """记录错误，支持链式调用"""
        self._items.append(FastqError(name, error_type, message))
        return self

    def __iter__(self) -> Iterator[FastqError]:
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)

    def __bool__(self) -> bool:
        return bool(self._items)

    def log_all(self, level: str = "error") -> None:
        """输出所有错误日志"""
        for err in self._items:
            err.log(level)

    def clear(self) -> None:
        """清空错误记录"""
        self._items.clear()


# ============================================================================
# 路径解析器 (PathResolver)
# ============================================================================


@dataclass
class PathResolver:
    """
    路径解析器 - 处理路径相关的所有操作

    遵循单一职责原则：只负责路径的加载、验证和匹配
    """

    include_paths: frozenset[Path] = field(default_factory=frozenset)
    exclude_paths: frozenset[Path] = field(default_factory=frozenset)

    def __post_init__(self) -> None:
        self._validate_no_overlap()

    def _validate_no_overlap(self) -> None:
        """验证包含和排除路径没有重叠"""
        if overlap := self.include_paths & self.exclude_paths:
            raise PathOverlapError(f"路径重叠: {', '.join(map(str, overlap))}")

    @classmethod
    def from_files(
        cls,
        include_file: Path | None = None,
        exclude_file: Path | None = None,
    ) -> Self:
        """从文件创建路径解析器"""
        return cls(
            include_paths=cls._load_path_set(include_file),
            exclude_paths=cls._load_path_set(exclude_file),
        )

    @staticmethod
    def _load_path_set(file_path: Path | None) -> frozenset[Path]:
        """从文件加载路径集合"""
        if file_path is None:
            return frozenset()

        if not file_path.exists():
            raise FileNotFoundError(f"路径列表文件不存在: {file_path}")

        paths = set()
        for line in filter(None, map(str.strip, file_path.read_text().splitlines())):
            resolved = Path(line).resolve()
            if not resolved.exists():
                raise FileNotFoundError(f"列表中的路径不存在: {resolved}")
            paths.add(resolved)

        return frozenset(paths)

    def is_excluded(self, path: Path) -> bool:
        """检查路径是否被排除"""
        resolved = path.resolve()
        return resolved in self.exclude_paths or any(resolved.is_relative_to(p) for p in self.exclude_paths)

    def find_target_dirs(
        self,
        base_dir: Path,
        target_names: Iterable[str],
        config: Config = DEFAULT_CONFIG,
    ) -> list[Path]:
        """查找目标目录"""
        target_set = set(target_names)
        targets = []

        # 从base_dir查找
        for date_dir in base_dir.glob(config.date_dir_pattern):
            if not date_dir.is_dir():
                continue
            for child in date_dir.iterdir():
                if child.name in target_set and not self.is_excluded(child):
                    targets.append(child.resolve())

        # 添加include路径中的目录
        for inc_path in self.include_paths:
            targets.extend(c.resolve() for c in inc_path.glob("*") if c.is_dir() and c.resolve() not in targets)

        return targets


# ============================================================================
# FASTQ扫描器 (FastqScanner)
# ============================================================================


class FastqScanner:
    """
    FASTQ文件扫描器 - 扫描目录并解析FASTQ文件信息

    遵循单一职责原则：只负责FASTQ文件的发现和解析
    """

    def __init__(self, errors: ErrorRecorder, config: Config = DEFAULT_CONFIG):
        self.errors = errors
        self.config = config

    def scan_directory(self, sample_path: Path) -> list[FastqInfo]:
        """扫描样品目录中的FASTQ文件"""
        if not sample_path.exists():
            logger.warning(f"样品路径不存在: {sample_path}")
            return []

        fastqs = self._glob_fastqs(sample_path)
        lib_id = self._extract_lib_id(sample_path)

        if not fastqs:
            logger.warning(f"在 {sample_path} 中未找到FASTQ文件")
            self.errors.record(lib_id, FastqErrorType.INCOMPLETE, f"未找到FASTQ文件 {sample_path}")
            return []

        results = []
        for fq in fastqs:
            if read_type := self._determine_read_type(fq.name):
                results.append(FastqInfo(lib_id, read_type, str(fq.absolute())))
            else:
                logger.warning(f"无法识别的FASTQ文件: {fq.name}")
                self.errors.record(lib_id, FastqErrorType.FORMAT, f"无法识别: {fq.name}")

        self._check_r1r2_pairing(fastqs, lib_id)
        return results

    def _glob_fastqs(self, path: Path) -> list[Path]:
        """获取目录下所有FASTQ文件"""
        return list(chain.from_iterable(path.glob(ext) for ext in self.config.fastq_extensions))

    def _extract_lib_id(self, lib_path: Path) -> str:
        """从路径提取文库ID"""
        parts = lib_path.name.split("-")
        lib_id = parts[-1]
        # 如果是数字或者单个字母（不区分大小写），则与前一部分组合
        if lib_id.isdigit() or (len(lib_id) == 1 and lib_id.isalpha()):
            lib_id = "-".join(parts[-2:])
        return lib_id

    def _determine_read_type(self, filename: str) -> str | None:
        """根据文件名确定读取类型(R1/R2)"""
        return next(
            (rt for pattern, rt in self.config.read_type_patterns.items() if filename.endswith(pattern)),
            None,
        )

    def _check_r1r2_pairing(self, fastqs: list[Path], lib_id: str) -> None:
        """检查R1/R2配对"""
        r1 = sum(1 for f in fastqs if "_R1.fastq.gz" in f.name)
        r2 = sum(1 for f in fastqs if "_R2.fastq.gz" in f.name)
        if r1 != r2:
            self.errors.record(lib_id, FastqErrorType.FORMAT, f"R1/R2不配对: {r1} != {r2}")


# ============================================================================
# 配置文件构建器 (ConfigBuilder)
# ============================================================================


class ConfigBuilder:
    """
    配置文件构建器 - 构建和管理libid-fastq配置

    遵循单一职责原则：只负责配置文件的读取、构建和保存
    """

    def __init__(
        self,
        scanner: FastqScanner,
        errors: ErrorRecorder,
        config: Config = DEFAULT_CONFIG,
    ):
        self.scanner = scanner
        self.errors = errors
        self.config = config
        self._line_tracker: list[LineTrack] = []
        self._duplicated_df = pd.DataFrame()

    @property
    def duplicated_data(self) -> pd.DataFrame:
        """获取重复数据"""
        return self._duplicated_df

    def build_or_load(self, fq_line_dir: Path, force_rebuild: bool = False) -> pd.DataFrame:
        """读取或构建配置文件"""
        config_file = fq_line_dir / self.config.config_filename
        self._track_directory(fq_line_dir)

        # 尝试读取现有配置
        if not force_rebuild and config_file.exists():
            if (df := self._try_load_config(config_file)) is not None:
                return df

        # 构建新配置
        return self._build_and_save_config(fq_line_dir, config_file)

    def _try_load_config(self, config_file: Path) -> pd.DataFrame | None:
        """尝试加载现有配置"""
        try:
            logger.info(f"使用已有配置: {config_file}")
            df = pd.read_table(config_file)
            if self.config.required_config_cols <= set(df.columns):
                self._validate_lib_r1r2(df, config_file.parent)
                return df
            logger.warning(f"配置缺少必要列，重建: {config_file}")
        except Exception as e:
            logger.warning(f"读取配置失败，重建: {e}")
        return None

    def _build_and_save_config(self, fq_line_dir: Path, config_file: Path) -> pd.DataFrame:
        """构建并保存新配置"""
        libid_map = self._build_libid_map(fq_line_dir)
        self._validate_lib_r1r2(libid_map, fq_line_dir)

        if not libid_map.empty:
            libid_map["dir_name"] = fq_line_dir.name
            try:
                libid_map.to_csv(config_file, sep="\t", index=False)
                logger.success(f"配置已保存: {config_file}")
            except Exception as e:
                logger.error(f"保存配置失败: {e}")

        return libid_map

    def _build_libid_map(self, fastq_path: Path) -> pd.DataFrame:
        """构建library ID到FASTQ文件的映射"""
        if not fastq_path.exists():
            logger.error(f"FASTQ路径不存在: {fastq_path}")
            return pd.DataFrame()

        sample_dirs = list(fastq_path.glob(self.config.sample_dir_pattern))
        if not sample_dirs:
            raise ValueError(f"在 {fastq_path} 中未找到{self.config.sample_dir_pattern}目录")

        records = [
            info._asdict()
            for path in tqdm(sample_dirs, desc=f"处理 {fastq_path.name}")
            for info in self.scanner.scan_directory(path)
        ]
        return pd.DataFrame(records)

    def _track_directory(self, fastq_path: Path) -> None:
        """跟踪目录结构"""
        self._line_tracker.extend(
            LineTrack(fastq_path.name, p.name, str(fastq_path.absolute()))
            for p in fastq_path.glob(self.config.sample_dir_pattern)
        )

    def _validate_lib_r1r2(self, df: pd.DataFrame, fq_line_dir: Path) -> None:
        """验证文库R1/R2配对"""
        if df.empty or "libid" not in df.columns or "read_type" not in df.columns:
            return

        counts = df[["libid", "read_type"]].value_counts().unstack(fill_value=0)
        if "R1" not in counts.columns or "R2" not in counts.columns:
            return

        for libid in counts[counts["R1"] != counts["R2"]].index:
            self.errors.record(
                str(libid),
                FastqErrorType.INCOMPLETE,
                f"{fq_line_dir} - {libid} R1/R2不配对",
            )

    def check_duplicates(self) -> None:
        """检查并记录重复数据"""
        if not self._line_tracker:
            return

        df = pd.DataFrame(self._line_tracker)
        dup_mask = df.duplicated(subset=["line", "lib_dir"], keep=False)  # type: ignore[arg-type]
        dup_df = df[dup_mask]

        if dup_df.empty:
            return

        for line, group in dup_df.groupby("line"):
            group_df = cast(pd.DataFrame, group)
            if len(group_df) <= 1:
                continue
            lib_dirs = cast(pd.Series, group_df["lib_dir"])
            line_paths = cast(pd.Series, group_df["line_path"])
            dirs = lib_dirs.unique()
            dir_names = ",".join(dirs[:3]) + (f" ...共{len(dirs)}个" if len(dirs) > 3 else "")
            paths = " | ".join(line_paths.unique())
            self.errors.record(
                str(line),
                FastqErrorType.DUPLICATED,
                f"<cyan>目录</cyan> {paths} <cyan>重复数据:</cyan> <w>{dir_names}</w>",
            )

        self._duplicated_df = (
            dup_df.groupby(["line", "lib_dir"])["line_path"]
            .agg(lambda s: " | ".join(pd.Series(s).unique()))
            .reset_index()
        )


# ============================================================================
# 样本验证器 (SampleValidator)
# ============================================================================


class SampleValidator:
    """
    样本验证器 - 验证样本信息的完整性和正确性

    遵循单一职责原则：只负责样本数据的验证
    """

    def __init__(
        self,
        errors: ErrorRecorder,
        warnings: ErrorRecorder,
        config: Config = DEFAULT_CONFIG,
    ):
        self.errors = errors
        self.warnings = warnings
        self.config = config

    def validate_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """验证必需列存在"""
        if missing := self.config.required_sample_cols - set(df.columns):
            raise ValidationError(f"缺少必要列: {missing}")
        return df

    def remove_duplicates(self, df: pd.DataFrame) -> pd.DataFrame:
        """去除重复行"""
        original = len(df)
        df = df.drop_duplicates()
        if (removed := original - len(df)) > 0:
            logger.warning(f"去除 {removed} 行重复")
        return df

    def check_nulls(self, df: pd.DataFrame) -> None:
        """检查空值"""
        nulls = df[list(self.config.required_sample_cols)].isnull().sum()
        if nulls.any():
            logger.warning(f"空值: {nulls[nulls > 0].to_dict()}")

    def check_sample_mapping(self, df: pd.DataFrame, threshold: float) -> None:
        """检查样本映射关系"""
        logger.info("检查样本映射...")
        self._check_duplicated_mapping(df)
        self._check_low_data(df, threshold)

    def _check_duplicated_mapping(self, df: pd.DataFrame) -> None:
        """检查重复映射"""
        dup_mask = df.duplicated(subset=["libid", "sample_id", "dir_name"])  # type: ignore[arg-type]
        dup_rows = df.loc[dup_mask, ["libid", "sample_id", "dir_name"]].to_dict("records")
        for row in dup_rows:
            self.warnings.record(
                str(row["libid"]),
                FastqErrorType.DUPLICATED,
                f"重复: {row['libid']} | {row['sample_id']} | {row['dir_name']}",
            )

    def _check_low_data(self, df: pd.DataFrame, threshold: float) -> None:
        """检查低数据量"""
        low_df = df[df["data_size"] < threshold]
        if not low_df.empty:
            logger.warning(f"{len(low_df)}个样本数据量 < {threshold}G")
            for row in low_df[["libid", "sample_id", "dir_name"]].to_dict("records"):
                self.errors.record(
                    str(row["libid"]),
                    FastqErrorType.INCOMPLETE,
                    f"数据量过低: {row['libid']} | {row['sample_id']} | {row['dir_name']}",
                )

    def validate(self, df: pd.DataFrame, threshold: float) -> pd.DataFrame:
        """完整验证流程"""
        df = self.validate_columns(df)
        df = self.remove_duplicates(df)
        self.check_nulls(df)
        self.check_sample_mapping(df, threshold)
        return df


# ============================================================================
# 脚本生成器 (ScriptGenerator)
# ============================================================================


class ScriptGenerator:
    """
    脚本生成器 - 生成合并/链接脚本

    遵循单一职责原则：只负责脚本文件的生成
    """

    def __init__(self, config: Config = DEFAULT_CONFIG):
        self.config = config

    def generate(
        self,
        fq_df: pd.DataFrame,
        output_dir: Path,
        errors: ErrorRecorder,
        mode: DataMode = DataMode.LINK,
    ) -> int:
        """生成脚本并返回生成数量"""
        scripts_dir = output_dir / "scripts"
        scripts_dir.mkdir(exist_ok=True, parents=True)

        script_count = 0
        for (sample_id, read_type), group in fq_df.groupby(["sample_id", "read_type"]):  # type: ignore[misc]
            if group["path"].isna().any():
                logger.warning(f"缺失路径: {sample_id}-{read_type}")
                errors.record(
                    str(sample_id),
                    FastqErrorType.INCOMPLETE,
                    f"缺失: {sample_id}-{read_type}",
                )
                continue

            script_file = self._create_script(
                scripts_dir,
                output_dir,
                sample_id,
                read_type,  # type: ignore[arg-type]
                sorted(group["path"].tolist()),
                mode,
            )
            if script_file:
                script_count += 1

        if script_count:
            logger.info(f"生成 {script_count} 个脚本")

        return script_count

    def _create_script(
        self,
        scripts_dir: Path,
        output_dir: Path,
        sample_id: str,
        read_type: str,
        fq_list: list[str],
        mode: DataMode,
    ) -> Path | None:
        """创建单个脚本文件"""
        out_fq = output_dir.absolute() / f"{sample_id}.{read_type}.fq.gz"
        cmd = self._make_command(fq_list, str(out_fq), mode)

        script_file = scripts_dir / f"{self.config.script_prefix}-{sample_id}-{read_type}.sh"
        script_file.write_text(f"#!/bin/bash\nset -euo pipefail\n{cmd}\n")
        script_file.chmod(0o755)
        return script_file

    @staticmethod
    def _make_command(fq_list: list[str], output: str, mode: DataMode) -> str:
        """生成合并或链接命令"""
        if len(fq_list) == 1:
            cmd = "ln -s" if mode == DataMode.LINK else "cp"
            return f"{cmd} {fq_list[0]} {output}"
        return f"cat {' '.join(fq_list)} > {output}"


# ============================================================================
# 脚本执行器 (ScriptRunner)
# ============================================================================


class ScriptRunner:
    """
    脚本执行器 - 并行执行脚本

    遵循单一职责原则：只负责脚本的执行
    """

    def __init__(self, config: Config = DEFAULT_CONFIG):
        self.config = config

    def run_parallel(self, scripts_dir: Path, max_workers: int = 8) -> dict[str, int]:
        """并行运行脚本"""
        pattern = f"{self.config.script_prefix}-*.sh"
        scripts = list(scripts_dir.glob(pattern))

        if not scripts:
            logger.warning(f"在 {scripts_dir} 中未找到脚本")
            return {"success": 0, "failed": 0}

        logger.info(f"并行执行 {len(scripts)} 个脚本")
        results = {"success": 0, "failed": 0}

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(self._run_script, s): s for s in scripts}
            for future in tqdm(as_completed(futures), total=len(futures), desc="执行脚本"):
                ok, _ = future.result()
                results["success" if ok else "failed"] += 1

        logger.info(f"结果: 成功 {results['success']}, 失败 {results['failed']}")
        return results

    @staticmethod
    def _run_script(script_path: Path) -> tuple[bool, str]:
        """运行单个脚本"""
        try:
            subprocess.run(["bash", str(script_path)], check=True, capture_output=True, text=True)
            return True, f"成功: {script_path.name}"
        except subprocess.CalledProcessError as e:
            msg = f"失败: {script_path.name} - {e.stderr}"
            logger.error(msg)
            return False, msg


# ============================================================================
# 统计日志器 (StatisticsLogger)
# ============================================================================


class StatisticsLogger:
    """统计信息日志器"""

    @staticmethod
    def log(df: pd.DataFrame) -> None:
        """记录统计信息"""
        if df.empty:
            logger.error("数据为空")
            return

        samples, libs = df["sample_id"].nunique(), df["libid"].nunique()
        logger.info(f"总计: {samples} 样品, {libs} 文库")

        miss = df[df["path"].isna()]
        if miss.empty:
            logger.success("所有数据已找到")
        else:
            logger.error(f"缺失: {miss['sample_id'].nunique()} 样品, {miss['libid'].nunique()} 文库")
            logger.error(f"示例: {list(miss['sample_id'].unique()[:10])}")


# ============================================================================
# 处理流水线 (Pipeline)
# ============================================================================


@dataclass
class Pipeline:
    """
    处理流水线 - 协调各组件完成整个处理流程

    遵循单一职责原则：只负责流程的编排和协调
    """

    sample_info: Path
    base_dir: Path
    output_dir: Path
    check_file: Path
    threads: int
    force_rebuild: bool
    threshold: float
    exclude: Path | None
    include: Path | None
    mode: DataMode
    run_script: bool
    config: Config = field(default_factory=Config)

    def execute(self) -> None:
        """执行处理流程"""
        # 1. 验证输入
        if not self.sample_info.exists():
            logger.error(f"样品文件不存在: {self.sample_info}")
            raise typer.Exit(1)

        # 2. 读取样品信息
        sample_df = self._load_sample_info()

        # 3. 初始化组件
        errors, warnings = ErrorRecorder(), ErrorRecorder()
        path_resolver = PathResolver.from_files(self.include, self.exclude)
        scanner = FastqScanner(errors, self.config)
        config_builder = ConfigBuilder(scanner, errors, self.config)
        validator = SampleValidator(errors, warnings, self.config)
        script_generator = ScriptGenerator(self.config)
        script_runner = ScriptRunner(self.config)

        # 配置日志文件
        logger.add(self.config.log_file, rotation="10 MB", retention="1 week", level="INFO")

        # 4. 验证样品信息
        sample_df = validator.validate(sample_df, self.threshold)
        sample_libs = sample_df["dir_name"].unique()
        logger.info(f"处理 {len(sample_libs)} 个目录")

        # 5. 去重
        sample_df = sample_df.drop_duplicates(subset=["sample_id", "dir_name", "libid"])  # type: ignore[call-overload]

        # 6. 查找目标目录
        target_dirs = path_resolver.find_target_dirs(self.base_dir, sample_libs, self.config)
        if not target_dirs:
            logger.error("未找到匹配目录")
            raise typer.Exit(1)

        logger.info(f"找到 {len(target_dirs)} 个匹配目录")

        # 7. 加载配置
        libid_map = self._load_all_configs(target_dirs, config_builder)
        if libid_map.empty:
            logger.error("未找到配置")
            raise typer.Exit(1)

        # 8. 合并数据
        logger.info("合并数据")
        merged = sample_df.merge(libid_map, how="left")
        StatisticsLogger.log(merged)

        # 9. 保存检查结果
        self._save_check_file(merged)

        # 10. 记录缺失数据
        self._record_missing_data(merged, errors)

        # 11. 生成输出
        self.output_dir.mkdir(exist_ok=True, parents=True)
        logger.info(f"输出到: {self.output_dir}")

        script_count = 0
        if self.run_script:
            script_count = script_generator.generate(merged, self.output_dir, errors, self.mode)

        # 12. 输出日志
        errors.log_all("error")
        warnings.log_all("warning")

        # 13. 保存重复数据
        self._save_duplicated_data(config_builder)

        # 14. 执行脚本
        if self.run_script:
            if errors:
                logger.warning("发现错误，跳过执行")
            elif script_count == 0:
                logger.warning("无脚本生成")
            else:
                results = script_runner.run_parallel(self.output_dir / "scripts", self.threads)
                logger.info(f"执行结果: {results}")
        else:
            if not errors and not warnings:
                logger.success("检查完成：无问题")

        logger.success(f"完成！输出: {self.output_dir}")

    def _load_sample_info(self) -> pd.DataFrame:
        """加载样品信息"""
        logger.info(f"读取样品信息: {self.sample_info}")
        try:
            col_map = self.config.sample_info_columns
            return pd.read_table(
                self.sample_info,
                header=None,
                names=list(col_map.keys()),
                usecols=list(col_map.values()),
            )
        except Exception as e:
            logger.error(f"读取失败: {e}")
            raise typer.Exit(1)

    def _load_all_configs(self, target_dirs: list[Path], builder: ConfigBuilder) -> pd.DataFrame:
        """加载所有配置"""
        logger.info("加载FASTQ配置")
        configs = []

        # 使用线程池并行加载配置
        with ThreadPoolExecutor(max_workers=self.config.max_threads) as executor:
            future_to_path = {
                executor.submit(builder.build_or_load, path, self.force_rebuild): path for path in target_dirs
            }

            for future in tqdm(as_completed(future_to_path), total=len(target_dirs), desc="加载配置"):
                path = future_to_path[future]
                try:
                    if not (cfg := future.result()).empty:
                        configs.append(cfg)
                except Exception as e:
                    logger.error(f"处理 {path} 出错: {e}")

        builder.check_duplicates()

        if not configs:
            return pd.DataFrame()

        result = pd.concat(configs, ignore_index=True)
        logger.success(f"成功加载 {len(result)} 条记录")
        return result

    def _save_check_file(self, df: pd.DataFrame) -> None:
        """保存检查结果"""
        try:
            df.to_csv(self.check_file, sep="\t", index=False)
            logger.success(f"检查结果: {self.check_file}")
        except Exception as e:
            logger.error(f"保存失败: {e}")

    def _record_missing_data(self, df: pd.DataFrame, errors: ErrorRecorder) -> None:
        """记录缺失数据"""
        for row in df[df["path"].isna()].itertuples():
            errors.record(
                str(row.libid),
                FastqErrorType.INCOMPLETE,
                f"{row.libid}-{row.sample_id}-{row.dir_name} 无数据",
            )

    def _save_duplicated_data(self, builder: ConfigBuilder) -> None:
        """保存重复数据"""
        if not builder.duplicated_data.empty:
            dup_file = self.output_dir / "duplicated_data.tsv"
            try:
                builder.duplicated_data.to_csv(dup_file, sep="\t", index=False)
                logger.success(f"重复数据: {dup_file}")
            except Exception as e:
                logger.error(f"保存重复数据失败: {e}")


# ============================================================================
# CLI 命令行接口
# ============================================================================


app = typer.Typer(help="FASTQ文件处理和合并工具", no_args_is_help=True)

# 参数类型定义
SampleInfoArg = Annotated[Path, typer.Argument(help="样品信息TSV文件")]
BaseDirOpt = Annotated[Path, typer.Option("--base-dir", "-b", help="FASTQ数据基础目录")]
OutputDirOpt = Annotated[Path, typer.Option("--output", "-o", help="输出目录")]
CheckFileOpt = Annotated[Path, typer.Option("--check-file", "-c", help="检查结果文件")]
ThreadsOpt = Annotated[int, typer.Option("--threads", "-t", min=1, max=32, help="并行线程数")]
ForceOpt = Annotated[bool, typer.Option("--force", "-f", help="强制重建配置")]
ThresholdOpt = Annotated[float, typer.Option("--threshold", help="数据量阈值(G)")]
ExcludeOpt = Annotated[Path | None, typer.Option("--exclude", "-e", help="排除路径列表文件")]
IncludeOpt = Annotated[Path | None, typer.Option("--include", "-i", help="包含路径列表文件")]
ModeOpt = Annotated[DataMode, typer.Option("--mode", "-m", help="处理模式")]


def _run_pipeline(pipeline: Pipeline) -> None:
    """执行流水线并处理异常"""
    try:
        pipeline.execute()
    except KeyboardInterrupt:
        logger.info("用户中断")
        raise typer.Exit(130)
    except typer.Exit:
        raise
    except Exception as e:
        logger.error(f"错误: {e}")
        raise typer.Exit(1)


@app.command()
def run(
    sample_info: SampleInfoArg,
    base_dir: BaseDirOpt = DEFAULT_CONFIG.base_dir,
    output_dir: OutputDirOpt = Path("."),
    check_file: CheckFileOpt = Path("check_file.tsv"),
    threads: ThreadsOpt = DEFAULT_CONFIG.default_threads,
    force: ForceOpt = False,
    threshold: ThresholdOpt = DEFAULT_CONFIG.default_threshold,
    exclude: ExcludeOpt = None,
    include: IncludeOpt = None,
    mode: ModeOpt = DataMode.LINK,
) -> None:
    """处理FASTQ文件并执行合并脚本"""
    _run_pipeline(
        Pipeline(
            sample_info,
            base_dir,
            output_dir,
            check_file,
            threads,
            force,
            threshold,
            exclude,
            include,
            mode,
            run_script=True,
        )
    )


@app.command()
def validate(
    sample_info: SampleInfoArg,
    base_dir: BaseDirOpt = DEFAULT_CONFIG.base_dir,
    output_dir: OutputDirOpt = Path("."),
    check_file: CheckFileOpt = Path("check_file.tsv"),
    threads: ThreadsOpt = DEFAULT_CONFIG.default_threads,
    force: ForceOpt = False,
    threshold: ThresholdOpt = DEFAULT_CONFIG.default_threshold,
    exclude: ExcludeOpt = None,
    include: IncludeOpt = None,
) -> None:
    """仅验证FASTQ文件，不执行合并脚本"""
    _run_pipeline(
        Pipeline(
            sample_info,
            base_dir,
            output_dir,
            check_file,
            threads,
            force,
            threshold,
            exclude,
            include,
            DataMode.LINK,
            run_script=False,
        )
    )


if __name__ == "__main__":
    app()
