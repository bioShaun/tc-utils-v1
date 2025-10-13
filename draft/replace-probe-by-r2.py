from pathlib import Path
from typing import Union

import pandas as pd
import typer


def kb_distance(a: int, b: int) -> float:
    """计算两个位置的千碱基距离，返回浮点数或整数"""
    distance = abs((b - a) / 1000)
    if distance < 1:
        return round(distance, 1)
    return round(distance, 0)


def main(
    ori_file: Path,
    r2_file: Path,
    patch_file: Path,
    out_file: Path,
) -> None:
    """
    Replace probes in origin file with probes from R2 file,
    by selecting the probe with the highest R2 for each target.

    Parameters
    ----------
    ori_file : Path
        Original probe file
    r2_file : Path
        R2 probe file
    patch_file : Path
        Patch file with target_id column
    out_file : Path
        Output file path

    Returns
    -------
    None
    """
    ori_df = pd.read_table(ori_file)
    r2_df = pd.read_table(r2_file, sep=r"\s+")
    patch_df = pd.read_table(patch_file)

    # 距离与列重命名
    r2_df["distance_to_origin_kb"] = r2_df.apply(
        lambda x: kb_distance(x["BP_A"], x["BP_B"]), axis=1
    )
    r2_df = r2_df.rename(
        columns={
            "CHR_A": "origin_chrom",
            "BP_A": "origin_pos",
            "SNP_A": "id_a",
            "CHR_B": "chrom",
            "BP_B": "pos",
            "SNP_B": "target_id",
        }
    )

    # 每个 target 选择最高 R²
    best_r2_df = r2_df.dropna(subset=["target_id"])
    best_r2_df = best_r2_df.loc[best_r2_df.groupby("target_id")["R2"].idxmax()]

    # 合并 patch
    add_r2_patch_df = patch_df.merge(
        best_r2_df.drop(["origin_chrom", "id_a", "chrom", "pos"], axis=1),
        on="target_id",
    )

    # 合并原始 probe level
    ori_level_df = ori_df[["chrom", "pos", "probe_level"]].drop_duplicates()
    ori_level_df.columns = ["chrom", "origin_pos", "origin_probe_level"]
    add_ori_level_df = add_r2_patch_df.merge(ori_level_df, on=["chrom", "origin_pos"])

    # 过滤：patch level < origin level
    filt_df = add_ori_level_df[
        add_ori_level_df["probe_level"] < add_ori_level_df["origin_probe_level"]
    ].copy()

    # R²分箱
    filt_df["R2_bin"] = (filt_df["R2"] / 0.05).astype(int)

    # 选出每个 origin 最优记录
    best_patch_df = filt_df.sort_values(
        ["R2_bin", "distance_to_origin_kb"], ascending=[False, True]
    ).drop_duplicates(subset=["chrom", "origin_pos"])

    # 标记并合并
    selected_patch_df = filt_df[
        filt_df["target_id"].isin(best_patch_df["target_id"])
    ].copy()
    selected_patch_df["ori_pos_id"] = (
        selected_patch_df["chrom"] + "_" + selected_patch_df["origin_pos"].astype(str)
    )
    if "pos_id" not in ori_df.columns:
        ori_df["pos_id"] = ori_df["chrom"].astype(str) + "_" + ori_df["pos"].astype(str)
    keep_ori_df = ori_df[~ori_df["pos_id"].isin(best_patch_df["ori_pos_id"])]

    merged_df = pd.concat(
        [
            keep_ori_df.drop(columns=["pos_id"]),
            selected_patch_df.drop(
                columns=["R2_bin", "ori_pos_id", "missing", "het", "maf"],
            ),
        ]
    )

    # 输出
    if out_file.suffix == ".xlsx":
        merged_df.to_excel(out_file, index=False)
    else:
        merged_df.to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)
