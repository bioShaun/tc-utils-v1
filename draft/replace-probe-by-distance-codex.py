#!/usr/bin/env python3
"""
根据距离替换 probe 的脚本。

以给定窗口大小（默认 10 kb）为步长，在最大范围（默认 100 kb）内
从候选表格中选择 maf 最高的 probe 进行替换。

需求点：
1. 在最大范围内找不到替换则直接报出所有无法替换的 probe 并退出；
2. 替换后的 probe 不允许重复；
3. 输出替换后的表格，并输出原始/替换位置的映射。
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import typer
from typing_extensions import Annotated


def _require_columns(df: pd.DataFrame, required: List[str], name: str) -> None:
    """Ensure required columns exist."""
    missing = [col for col in required if col not in df.columns]
    if missing:
        typer.echo(
            f"[{name}] 缺少必要列: {', '.join(missing)}",
            err=True,
        )
        raise typer.Exit(1)


def _prepare_df(df: pd.DataFrame, name: str) -> pd.DataFrame:
    """Cast key columns to expected types."""
    df = df.copy()
    df["chrom"] = df["chrom"].astype(str)
    try:
        df["pos"] = pd.to_numeric(df["pos"], errors="raise").astype(int)
    except Exception as exc:  # pragma: no cover - 类型转换错误时直接退出
        typer.echo(f"[{name}] pos 列无法转换为整数: {exc}", err=True)
        raise typer.Exit(1)

    # maf 只用于排序，无法转换时设置为 -1 以降低优先级
    df["maf"] = pd.to_numeric(df["maf"], errors="coerce").fillna(-1.0)

    df["pos_id"] = df["chrom"] + "_" + df["pos"].astype(str)
    return df


def _select_candidate(
    target: pd.Series,
    candidate_df: pd.DataFrame,
    used_pos_ids: set,
    window_bp: int,
    max_distance_bp: int,
) -> Optional[pd.Series]:
    """
    在候选表中为单个 probe 选择最佳替换。

    策略：按 window_bp 归类窗口，优先选择窗口索引最小且 maf 最大的候选；
    同一窗口内次要按距离、id 排序。
    """
    chrom_df = candidate_df[candidate_df["chrom"] == target["chrom"]]
    if chrom_df.empty:
        return None

    # 计算距离并限制最大范围
    chrom_df = chrom_df.assign(
        distance_bp=(chrom_df["pos"] - target["pos"]).abs()
    )
    chrom_df = chrom_df[chrom_df["distance_bp"] <= max_distance_bp]
    if chrom_df.empty:
        return None

    # 排除已使用的候选
    chrom_df = chrom_df[~chrom_df["pos_id"].isin(used_pos_ids)]
    if chrom_df.empty:
        return None

    chrom_df = chrom_df.assign(window_index=(chrom_df["distance_bp"] // window_bp))

    ranked = chrom_df.sort_values(
        ["window_index", "maf", "distance_bp", "id"],
        ascending=[True, False, True, True],
    )
    return ranked.iloc[0]


def _precheck_availability(
    replace_df: pd.DataFrame,
    candidate_df: pd.DataFrame,
    window_bp: int,
    max_distance_bp: int,
) -> Tuple[bool, pd.DataFrame]:
    """
    检查每个需要替换的 probe 是否能找到候选。

    使用与正式替换相同的去重逻辑（used_pos_ids），保证检查结果可靠。
    """
    used: set = set()
    missing_rows: List[pd.Series] = []

    for _, row in replace_df.iterrows():
        selected = _select_candidate(
            row,
            candidate_df,
            used,
            window_bp=window_bp,
            max_distance_bp=max_distance_bp,
        )
        if selected is None:
            missing_rows.append(row)
        else:
            used.add(selected["pos_id"])

    return len(missing_rows) == 0, pd.DataFrame(missing_rows)


def main(
    replace_file: Annotated[Path, typer.Argument(help="需要替换的 probe 表格")],
    candidate_file: Annotated[Path, typer.Argument(help="候选 probe 表格")],
    out_file: Annotated[Path, typer.Argument(help="输出替换结果的文件路径")],
    window_kb: Annotated[
        int,
        typer.Option(
            "--window-kb",
            "-w",
            help="搜索窗口大小（kb，默认 10kb）",
        ),
    ] = 10,
    max_distance_kb: Annotated[
        int,
        typer.Option(
            "--max-distance-kb",
            "-m",
            help="搜索最大范围（kb，默认 100kb）",
        ),
    ] = 100,
    mapping_file: Annotated[
        Optional[Path],
        typer.Option(
            "--mapping",
            "-p",
            help="替换映射输出文件（原始/新 probe 对应关系），默认与 out_file 同目录",
        ),
    ] = None,
) -> None:
    """
    根据距离替换 probe，窗口内优先 maf 高的候选。

    步骤：
    1) 检查所有需要替换的 probe 是否能在最大范围内找到候选，否则直接退出；
    2) 按窗口优先级选取候选，保证替换后的 probe 不重复；
    3) 输出替换后的表格及映射文件。
    """
    required_cols = ["chrom", "pos", "id", "target_id", "maf"]

    replace_df = pd.read_table(replace_file)
    candidate_df = pd.read_table(candidate_file)

    _require_columns(replace_df, required_cols, "replace")
    _require_columns(candidate_df, required_cols, "candidate")

    replace_df = _prepare_df(replace_df, "replace")
    candidate_df = _prepare_df(candidate_df, "candidate")

    # 删除候选表中与待替换表重复的条目（避免自替换）
    candidate_df = candidate_df[~candidate_df["pos_id"].isin(replace_df["pos_id"])]

    window_bp = window_kb * 1000
    max_distance_bp = max_distance_kb * 1000

    typer.echo(f"需要替换的 probe 数量: {len(replace_df)}")
    typer.echo(f"候选 probe 数量: {len(candidate_df)}")
    typer.echo(f"窗口大小: {window_kb} kb, 最大范围: {max_distance_kb} kb")

    ok, missing_df = _precheck_availability(
        replace_df, candidate_df, window_bp, max_distance_bp
    )
    if not ok:
        typer.echo(
            f"\n错误：有 {len(missing_df)} 个 probe 在 {max_distance_kb} kb 内找不到可用替换：",
            err=True,
        )
        cols_to_show = [col for col in ["chrom", "pos", "id", "target_id", "maf"] if col in missing_df.columns]
        typer.echo(missing_df[cols_to_show].to_string(index=False), err=True)

        no_replace_path = out_file.parent / f"{out_file.stem}_no_replacement.tsv"
        missing_df.drop(columns=["pos_id"], errors="ignore").to_csv(
            no_replace_path, sep="\t", index=False
        )
        typer.echo(f"无法替换的 probe 已输出到: {no_replace_path}", err=True)
        raise typer.Exit(1)

    # 正式替换
    used: set = set()
    replacements: List[pd.Series] = []
    mapping_records: List[Dict[str, object]] = []

    for _, row in replace_df.iterrows():
        selected = _select_candidate(
            row,
            candidate_df,
            used,
            window_bp=window_bp,
            max_distance_bp=max_distance_bp,
        )
        if selected is None:
            # 理论上不会发生（前置检查已通过），但仍做保护
            typer.echo(
                f"警告：{row['chrom']}:{row['pos']} 未找到可用替换，跳过。",
                err=True,
            )
            continue

        used.add(selected["pos_id"])

        selected = selected.copy()
        selected["origin_chrom"] = row["chrom"]
        selected["origin_pos"] = row["pos"]
        selected["origin_id"] = row["id"]
        selected["origin_target_id"] = row["target_id"]
        selected["distance_to_origin_bp"] = int(
            abs(int(selected["pos"]) - int(row["pos"]))
        )
        selected["window_index"] = int(selected["distance_to_origin_bp"] // window_bp)
        replacements.append(selected)

        mapping_records.append(
            {
                "origin_chrom": row["chrom"],
                "origin_pos": row["pos"],
                "origin_id": row["id"],
                "origin_target_id": row["target_id"],
                "replacement_chrom": selected["chrom"],
                "replacement_pos": selected["pos"],
                "replacement_id": selected["id"],
                "replacement_target_id": selected.get("target_id", None),
                "replacement_maf": selected["maf"],
                "distance_bp": selected["distance_to_origin_bp"],
                "window_index": selected["window_index"],
            }
        )

    if not replacements:
        typer.echo("未生成任何替换结果。", err=True)
        raise typer.Exit(1)

    replacement_df = pd.DataFrame(replacements)

    # 输出列：保留候选表原列，加上原始位点和距离信息
    base_cols = [col for col in candidate_df.columns if col != "pos_id"]
    extra_cols = [
        "origin_chrom",
        "origin_pos",
        "origin_id",
        "origin_target_id",
        "distance_to_origin_bp",
        "window_index",
    ]
    output_cols = [col for col in base_cols + extra_cols if col in replacement_df.columns]
    output_df = replacement_df[output_cols]

    if out_file.suffix == ".xlsx":
        output_df.to_excel(out_file, index=False)
    else:
        output_df.to_csv(out_file, sep="\t", index=False)
    typer.echo(f"替换结果已保存到: {out_file}")

    if mapping_file is None:
        mapping_file = out_file.parent / f"{out_file.stem}_mapping.tsv"

    mapping_df = pd.DataFrame(mapping_records)
    mapping_df.to_csv(mapping_file, sep="\t", index=False)
    typer.echo(f"替换映射已保存到: {mapping_file}")


if __name__ == "__main__":
    typer.run(main)
