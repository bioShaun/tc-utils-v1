#!/usr/bin/env python3
"""
根据距离替换probe的脚本。

以指定窗口大小（默认10kb）为步长，在最大范围（默认100kb）内
从候选probe中选择maf最高的probe进行替换。
"""

from pathlib import Path
from typing import Optional

import pandas as pd
import typer
from typing_extensions import Annotated


def find_best_replacement(
    origin_chrom: str,
    origin_pos: int,
    candidate_df: pd.DataFrame,
    used_candidates: set,
    window_kb: float,
    max_distance_kb: float,
) -> Optional[pd.Series]:
    """
    为一个需要替换的probe找到最佳替换候选。

    策略：以window_kb为步长，从近到远搜索，在每个窗口内选择maf最高的候选。
    确保选择的候选不在已使用的集合中。

    Parameters
    ----------
    origin_chrom : str
        原始probe的染色体
    origin_pos : int
        原始probe的位置
    candidate_df : pd.DataFrame
        候选probe表格
    used_candidates : set
        已经使用过的候选probe的pos_id集合
    window_kb : float
        窗口大小（kb）
    max_distance_kb : float
        最大搜索距离（kb）

    Returns
    -------
    Optional[pd.Series]
        最佳替换候选的行，如果找不到则返回None
    """
    # 筛选同一染色体的候选
    same_chrom_df = candidate_df[candidate_df["chrom"] == origin_chrom].copy()
    if same_chrom_df.empty:
        return None

    # 计算距离
    same_chrom_df["distance"] = abs(same_chrom_df["pos"] - origin_pos)
    same_chrom_df["distance_kb"] = same_chrom_df["distance"] / 1000

    # 只保留最大距离范围内的候选
    in_range_df = same_chrom_df[same_chrom_df["distance_kb"] <= max_distance_kb]
    if in_range_df.empty:
        return None

    # 排除已使用的候选
    in_range_df = in_range_df[~in_range_df["pos_id"].isin(used_candidates)]
    if in_range_df.empty:
        return None

    # 按窗口从近到远搜索
    window_size = window_kb * 1000  # 转换为bp
    current_start = 0

    while current_start < max_distance_kb * 1000:
        current_end = current_start + window_size
        window_df = in_range_df[
            (in_range_df["distance"] >= current_start)
            & (in_range_df["distance"] < current_end)
        ]

        if not window_df.empty:
            # 在当前窗口内选择maf最高的
            best_idx = window_df["maf"].idxmax()
            return in_range_df.loc[best_idx]

        current_start = current_end

    return None


def main(
    replace_file: Annotated[Path, typer.Argument(help="需要替换的probe表格文件")],
    candidate_file: Annotated[Path, typer.Argument(help="候选probe表格文件")],
    out_file: Annotated[Path, typer.Argument(help="输出文件路径")],
    window_kb: Annotated[
        float, typer.Option("--window", "-w", help="搜索窗口大小（kb）")
    ] = 10.0,
    max_distance_kb: Annotated[
        float, typer.Option("--max-distance", "-m", help="最大搜索距离（kb）")
    ] = 100.0,
    mapping_file: Annotated[
        Optional[Path],
        typer.Option(
            "--mapping", "-p", help="替换映射输出文件（原pos到新pos的对应关系）"
        ),
    ] = None,
) -> None:
    """
    根据距离替换probe。

    以window_kb为窗口，在max_distance_kb范围内从候选表格中
    选择maf最高的probe进行替换。

    替换前会检查是否所有需要替换的probe都能找到替换，
    如果有找不到替换的probe，会输出这些probe并停止。
    """
    # 读取数据
    replace_df = pd.read_table(replace_file)
    candidate_df = pd.read_table(candidate_file)

    # 添加pos_id用于标识
    replace_df["pos_id"] = (
        replace_df["chrom"].astype(str) + "_" + replace_df["pos"].astype(str)
    )
    candidate_df["pos_id"] = (
        candidate_df["chrom"].astype(str) + "_" + candidate_df["pos"].astype(str)
    )

    # 排除需要替换的probe本身（如果候选表格中包含）
    candidate_df = candidate_df[~candidate_df["pos_id"].isin(replace_df["pos_id"])]

    typer.echo(f"需要替换的probe数量: {len(replace_df)}")
    typer.echo(f"候选probe数量: {len(candidate_df)}")
    typer.echo(f"窗口大小: {window_kb} kb")
    typer.echo(f"最大搜索距离: {max_distance_kb} kb")

    # 第一轮：检查是否所有需要替换的probe都能找到替换
    typer.echo("\n检查是否所有probe都能找到替换...")
    cannot_replace = []
    used_candidates_check: set = set()

    for _, row in replace_df.iterrows():
        best = find_best_replacement(
            row["chrom"],
            row["pos"],
            candidate_df,
            used_candidates_check,
            window_kb,
            max_distance_kb,
        )
        if best is None:
            cannot_replace.append(row)
        else:
            used_candidates_check.add(best["pos_id"])

    if cannot_replace:
        typer.echo(
            f"\n错误: 有 {len(cannot_replace)} 个probe在{max_distance_kb}kb范围内找不到替换:",
            err=True,
        )
        cannot_replace_df = pd.DataFrame(cannot_replace)
        typer.echo(
            cannot_replace_df[["chrom", "pos", "id", "target_id", "maf"]].to_string(),
            err=True,
        )

        # 输出找不到替换的probe到文件
        no_replace_file = out_file.parent / f"{out_file.stem}_no_replacement.tsv"
        cannot_replace_df.drop(columns=["pos_id"]).to_csv(
            no_replace_file, sep="\t", index=False
        )
        typer.echo(f"\n找不到替换的probe已保存到: {no_replace_file}", err=True)
        raise typer.Exit(1)

    typer.echo("检查通过，所有probe都能找到替换。")

    # 第二轮：正式替换
    typer.echo("\n开始替换...")
    used_candidates: set = set()
    replacements = []
    mapping_records = []

    for _, row in replace_df.iterrows():
        best = find_best_replacement(
            row["chrom"],
            row["pos"],
            candidate_df,
            used_candidates,
            window_kb,
            max_distance_kb,
        )

        if best is not None:
            used_candidates.add(best["pos_id"])
            # 记录替换信息
            replacement_row = best.copy()
            replacement_row["origin_pos"] = row["pos"]
            replacement_row["origin_id"] = row["id"]
            replacement_row["distance_to_origin_kb"] = round(
                abs(best["pos"] - row["pos"]) / 1000, 2
            )
            replacements.append(replacement_row)

            mapping_records.append(
                {
                    "chrom": row["chrom"],
                    "origin_pos": row["pos"],
                    "origin_id": row["id"],
                    "new_pos": best["pos"],
                    "new_id": best["id"],
                    "distance_kb": round(abs(best["pos"] - row["pos"]) / 1000, 2),
                    "new_maf": best["maf"],
                }
            )

    # 创建替换后的DataFrame
    replacement_df = pd.DataFrame(replacements)

    # 整理输出列（保留原始表格的列，加上origin_pos和distance_to_origin_kb）
    original_cols = [col for col in candidate_df.columns if col != "pos_id"]
    output_cols = original_cols + ["origin_pos", "origin_id", "distance_to_origin_kb"]
    output_cols = [col for col in output_cols if col in replacement_df.columns]

    output_df = replacement_df[output_cols]

    # 输出替换后的probe
    if out_file.suffix == ".xlsx":
        output_df.to_excel(out_file, index=False)
    else:
        output_df.to_csv(out_file, sep="\t", index=False)

    typer.echo(f"\n替换完成，结果已保存到: {out_file}")
    typer.echo(f"成功替换的probe数量: {len(output_df)}")

    # 输出映射文件
    if mapping_file is None:
        mapping_file = out_file.parent / f"{out_file.stem}_mapping.tsv"

    mapping_df = pd.DataFrame(mapping_records)
    mapping_df.to_csv(mapping_file, sep="\t", index=False)
    typer.echo(f"替换映射已保存到: {mapping_file}")


if __name__ == "__main__":
    typer.run(main)
