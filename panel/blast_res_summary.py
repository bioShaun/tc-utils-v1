from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import Iterator, List, Optional

import pandas as pd
import typer
from tqdm import tqdm


@dataclass
class BlastConfig:
    """BLAST 分析配置类"""

    min_match_length: int = 24
    max_match_count: int = 1
    probe_length: int = 120
    output_all: bool = False
    show_max_length: bool = False
    min_identity: Optional[int] = None
    max_mismatch_count: Optional[int] = None
    max_gap_count: Optional[int] = None


def get_real_match_length(row: pd.Series, probe_length: int) -> int:
    """计算实际匹配长度

    Args:
        row: 包含匹配信息的 Series
        probe_length: 探针长度

    Returns:
        实际匹配长度
    """
    match_len = row["match_len"]
    mismatch = row["mismatch"]
    gapopen = row["gapopen"]

    return (
        probe_length - mismatch + gapopen
        if match_len > probe_length
        else match_len - mismatch - gapopen
    )


def get_blast_files(blast_dir: Path) -> List[Path]:
    """获取 BLAST 结果文件列表

    Args:
        blast_dir: BLAST 结果文件目录

    Returns:
        排序后的文件列表

    Raises:
        typer.Exit: 目录不存在或为空时退出
    """
    if not blast_dir.exists():
        typer.echo(f"目录不存在: {blast_dir}", err=True)
        raise typer.Exit(1)

    files = sorted(blast_dir.glob("*"))
    if not files:
        typer.echo(f"目录为空: {blast_dir}", err=True)
        raise typer.Exit(1)

    return files


def process_blast_file(
    blast_file: Path, config: BlastConfig, chunk_size: int = 10000
) -> pd.DataFrame:
    """处理单个 BLAST 文件

    Args:
        blast_file: BLAST 文件路径
        config: 分析配置
        chunk_size: 数据块大小

    Returns:
        处理后的数据框
    """
    # 使用 chunks 分块读取大文件
    chunks = pd.read_table(
        blast_file,
        usecols=[0, 2, 3, 4, 5],
        names=["id", "identity", "match_len", "mismatch", "gapopen"],
        chunksize=chunk_size,
    )

    blast_dfs = []
    for chunk in chunks:
        if not config.output_all:
            chunk = filter_blast_results(chunk, config)
        blast_dfs.append(chunk)

    return pd.concat(blast_dfs) if blast_dfs else pd.DataFrame()


def filter_blast_results(df: pd.DataFrame, config: BlastConfig) -> pd.DataFrame:
    """根据配置过滤 BLAST 结果

    Args:
        df: BLAST 结果数据框
        config: 分析配置

    Returns:
        过滤后的数据框
    """
    mask = df["match_len"] >= config.min_match_length

    if config.min_identity is not None:
        mask &= df["identity"] >= config.min_identity
    if config.max_mismatch_count is not None:
        mask &= df["mismatch"] <= config.max_mismatch_count
    if config.max_gap_count is not None:
        mask &= df["gapopen"] <= config.max_gap_count

    return df[mask]


def calculate_match_statistics(
    blast_df: pd.DataFrame, config: BlastConfig
) -> pd.DataFrame:
    """计算匹配统计信息

    Args:
        blast_df: BLAST 结果数据框
        config: 分析配置

    Returns:
        包含统计信息的数据框
    """
    if blast_df.empty:
        return pd.DataFrame(columns=["id", "blast_match", "second_max_length"])

    # 计算匹配次数
    id_count_df = blast_df["id"].value_counts()
    if not config.output_all:
        id_count_df = id_count_df[id_count_df <= config.max_match_count]
    id_count_df = id_count_df.reset_index()
    id_count_df.columns = ["id", "blast_match"]

    # 计算实际匹配长度
    my_get_real_match_length = partial(
        get_real_match_length, probe_length=config.probe_length
    )
    blast_df = blast_df.groupby("id").head(2).copy()
    blast_df["real_match_length"] = blast_df.apply(my_get_real_match_length, axis=1)

    # 计算次优匹配长度
    best_match_idx = blast_df.groupby("id")["real_match_length"].idxmax()
    other_match_df = blast_df[~blast_df.index.isin(best_match_idx)]
    second_max_length_df = (
        other_match_df.groupby("id")["real_match_length"]
        .max()
        .reset_index()
        .rename(columns={"real_match_length": "second_max_length"})
    )

    # 合并结果
    result_df = pd.merge(id_count_df, second_max_length_df, on="id", how="left")
    result_df["second_max_length"] = (
        result_df["second_max_length"].fillna(0).astype(int)
    )

    if config.show_max_length:
        max_length_df = (
            blast_df.groupby("id")["real_match_length"]
            .max()
            .reset_index()
            .rename(columns={"real_match_length": "max_length"})
        )
        result_df = pd.merge(result_df, max_length_df, on="id")

    return result_df


def main(
    blast_dir: Path = typer.Argument(..., help="BLAST 结果文件目录"),
    out_file: Path = typer.Argument(..., help="输出文件路径"),
    min_match_length: int = typer.Option(24, help="最小匹配长度"),
    max_match_count: int = typer.Option(1, help="最大匹配次数"),
    probe_length: int = typer.Option(120, help="探针长度"),
    output_all: bool = typer.Option(False, help="输出所有结果"),
    show_max_length: bool = typer.Option(False, help="显示最大长度"),
    min_identity: Optional[int] = typer.Option(None, help="最小一致性"),
    max_mismatch_count: Optional[int] = typer.Option(None, help="最大错配数"),
    max_gap_count: Optional[int] = typer.Option(None, help="最大间隔数"),
) -> None:
    """处理 BLAST 结果文件并生成统计报告"""
    try:
        # 初始化配置
        config = BlastConfig(
            min_match_length=min_match_length,
            max_match_count=max_match_count,
            probe_length=probe_length,
            output_all=output_all,
            show_max_length=show_max_length,
            min_identity=min_identity,
            max_mismatch_count=max_mismatch_count,
            max_gap_count=max_gap_count,
        )

        # 获取文件列表
        blast_files = get_blast_files(blast_dir)

        # 处理文件
        with tqdm(blast_files, desc="Processing BLAST files") as pbar:
            for n, blast_file in enumerate(pbar):
                pbar.set_postfix({"file": blast_file.name})

                # 处理单个文件
                blast_df = process_blast_file(blast_file, config)
                result_df = calculate_match_statistics(blast_df, config)

                # 写入结果
                mode = "w" if n == 0 else "a"
                header = n == 0
                result_df.to_csv(
                    out_file, sep="\t", index=False, mode=mode, header=header
                )

        typer.echo(f"结果已保存至: {out_file}")

    except Exception as e:
        typer.echo(f"处理过程中出现错误: {str(e)}", err=True)
        raise typer.Exit(1)


if __name__ == "__main__":
    typer.run(main)
