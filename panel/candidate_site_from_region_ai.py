import logging
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import pandas as pd
import typer
from rich.logging import RichHandler
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeRemainingColumn,
)

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, markup=True)],
)
logger = logging.getLogger(__name__)


@dataclass
class PriorityOrder:
    """
    PriorityOrder类用于定义优先级排序的规则。

    该类包含两个类属性：`columns`和`ascending`，分别指定了排序的列名和对应的排序顺序（升序或降序）。
    - `columns`: 包含排序依据的列名列表，默认为["sequence_score", "region_score", "maf"]。
    - `ascending`: 包含与`columns`对应的排序顺序列表，默认为[False, False, False]，表示所有列均按降序排序。
    """

    columns: List[str] = field(
        default_factory=lambda: ["sequence_score", "region_score", "maf"]
    )
    ascending: List[bool] = field(default_factory=lambda: [False, False, False])

    def __post_init__(self) -> None:
        """验证columns和ascending长度是否匹配"""
        if len(self.columns) != len(self.ascending):
            raise ValueError("columns和ascending的长度必须相同")


def validate_input_files(
    region_file: Path, candidate_site_file: Path, output_file: Path
) -> None:
    """验证输入输出文件的有效性"""
    if not region_file.exists():
        raise FileNotFoundError(f"Region file not found: {region_file}")
    if not candidate_site_file.exists():
        raise FileNotFoundError(f"Candidate site file not found: {candidate_site_file}")
    if output_file.exists():
        logger.warning(
            f"Output file {output_file} already exists and will be overwritten"
        )

    output_file.parent.mkdir(parents=True, exist_ok=True)


def process_region(
    region: pd.Series,
    sites_df: pd.DataFrame,
    region_size: int,
    candidate_per_region: int,
) -> Optional[pd.DataFrame]:
    """处理单个区域的数据"""
    try:
        start, end = int(region.start), int(region.end)

        region_sites = sites_df[
            (sites_df["chrom"] == region.chr)
            & (sites_df["pos"] > start)
            & (sites_df["pos"] <= end)
        ].copy()

        if region_sites.empty:
            return None

        region_sites["pos_range"] = pd.cut(
            region_sites["pos"], bins=range(start, end + 1, region_size), right=False
        )
        out_df = (
            region_sites.groupby("pos_range", observed=True)
            .head(candidate_per_region)
            .drop(columns=["pos_range"])
        )
        out_df["region"] = f"{region.chr}:{start}-{end}"

        return out_df
    except Exception as e:
        logger.error(f"Error processing region {region.chr}:{start}-{end}: {str(e)}")
        return None


def group_sites_by_region(df: pd.DataFrame) -> pd.DataFrame:
    id_region_df = (
        df.groupby("id")["region"].unique().map(lambda x: ",".join(x)).reset_index()
    )
    rm_dup_df = df.drop_duplicates(subset=["id"]).drop(columns=["region"])
    return id_region_df.merge(rm_dup_df, on="id")


def candidate_site_from_region(
    region_file: Path = typer.Argument(..., help="Path to the region file"),
    candidate_site_file: Path = typer.Argument(
        ..., help="Path to the candidate site file"
    ),
    output_file: Path = typer.Argument(..., help="Path to save the output"),
    max_region: int = typer.Option(
        100, "--max-region", "-m", help="Maximum number of regions to consider"
    ),
    candidate_per_region: int = typer.Option(
        1,
        "--candidate-per-region",
        "-c",
        help="Number of candidates to select per region",
    ),
    default_region_size: int = typer.Option(
        200, "--region-size", "-s", help="Default size of each region"
    ),
    threads: int = typer.Option(
        4, "--threads", "-t", help="Number of threads to use for processing"
    ),
) -> None:
    """
    Select candidate sites from regions based on specified criteria.

    This tool processes region and candidate site files to select optimal candidate sites
    based on various parameters and scoring criteria.
    """
    try:
        # 验证输入文件
        validate_input_files(region_file, candidate_site_file, output_file)

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TimeRemainingColumn(),
        ) as progress:
            # 读取数据文件
            file_task = progress.add_task("[cyan]Reading input files...", total=2)

            region_df = pd.read_csv(
                region_file,
                sep="\t",
                header=None,
                names=["chr", "start", "end"],
                usecols=[0, 1, 2],
                dtype={"chr": str, "start": int, "end": int},
            )
            progress.advance(file_task)

            sites_df = pd.read_csv(
                candidate_site_file,
                sep="\t",
                dtype={
                    "chrom": str,
                    "pos": int,
                    "sequence_score": float,
                    "region_score": float,
                    "maf": float,
                },
            )
            progress.advance(file_task)

            # 验证必需的列
            required_columns = ["chrom", "pos"] + PriorityOrder().columns
            missing_columns = [
                col for col in required_columns if col not in sites_df.columns
            ]
            if missing_columns:
                raise ValueError(
                    f"Missing required columns in candidate site file: {missing_columns}"
                )

            # 排序候选位点
            sort_task = progress.add_task("[yellow]Sorting candidate sites...", total=1)
            sites_df.sort_values(
                by=PriorityOrder().columns,
                ascending=PriorityOrder().ascending,
                inplace=True,
            )
            progress.advance(sort_task)

            # 处理每个区域
            total_regions = len(region_df)
            process_task = progress.add_task(
                "[green]Processing regions...", total=total_regions
            )
            select_dfs = []

            with ThreadPoolExecutor(max_workers=threads) as executor:
                futures = []
                for region in region_df.itertuples():
                    span = int(region.end) - int(region.start)
                    region_size = (
                        default_region_size
                        if default_region_size * max_region >= span
                        else span // max_region
                    )

                    futures.append(
                        executor.submit(
                            process_region,
                            region,
                            sites_df,
                            region_size,
                            candidate_per_region,
                        )
                    )

                for future in as_completed(futures):
                    result = future.result()
                    if result is not None:
                        select_dfs.append(result)
                    progress.advance(process_task)

            # 合并结果并保存
            save_task = progress.add_task("[blue]Saving results...", total=1)
            if not select_dfs:
                logger.warning("No candidate sites were selected!")
                return

            result_df = pd.concat(select_dfs, ignore_index=True)
            result_df = group_sites_by_region(result_df)
            result_df.to_csv(output_file, sep="\t", index=False)
            progress.advance(save_task)

            # 显示结果统计
            logger.info(f"Selected {len(result_df)} candidate sites")
            logger.info(f"Results saved to {output_file}")

    except Exception as e:
        logger.error(f"An error occurred: {str(e)}")
        sys.exit(1)


def main():
    """Entry point for the application."""
    typer.run(candidate_site_from_region)


if __name__ == "__main__":
    main()
