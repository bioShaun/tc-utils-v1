#!/usr/bin/env python
# -*- coding: utf-8 -*-"""
"""
从gffcompare的.tmap文件中提取基因ID映射 (v2, 使用pandas优化)

规则:
1. class code 为 =, c, k, m, j 认为是同一基因，保留所有配对
2. 可靠度排序: = > c ≈ k > m > j
"""

import argparse
import sys
import pandas as pd


def get_class_priority(class_code_series):
    """返回class code的优先级，数值越小优先级越高"""
    priority_map = {"=": 1, "c": 2, "k": 2, "m": 3, "j": 4}
    return class_code_series.map(priority_map).fillna(999).astype(int)


def parse_tmap_pandas(tmap_file, valid_classes={"=", "c", "k", "m", "j"}):
    """
    使用pandas解析tmap文件，提取基因ID映射

    tmap文件格式(tab分隔):
    ref_gene_id  class_code  qry_gene_id  ...
    """
    # 尝试读取header来确定列名
    with open(tmap_file, "r") as f:
        header = f.readline().strip().split("\t")

    # 检查关键列是否存在
    if "ref_gene_id" in header and "class_code" in header and "qry_gene_id" in header:
        use_cols = ["ref_gene_id", "class_code", "qry_gene_id"]
        df = pd.read_csv(
            tmap_file, sep="\t", usecols=use_cols, na_values=["-"], comment="#"
        )
    else:
        # 如果header不规范，则按位置读取，并假设没有header
        # 标准tmap格式: 列1是ref_gene_id, 列3是class_code, 列5是qry_gene_id
        try:
            df = pd.read_csv(
                tmap_file, sep="\t", header=None, na_values=["-"], comment="#"
            )
            # 提取需要的列
            df = df[[0, 2, 4]]
            df.columns = ["ref_gene_id", "class_code", "qry_gene_id"]
        except (IndexError, KeyError):
            # 如果列数不足，返回空DataFrame
            return pd.DataFrame(columns=["qry_gene_id", "ref_gene_id", "class_code"])

    # 删除包含空值的行
    df.dropna(inplace=True)

    # 过滤无效的class code
    df = df[df["class_code"].isin(valid_classes)]

    if df.empty:
        return pd.DataFrame(columns=["qry_gene_id", "ref_gene_id", "class_code"])

    # 计算优先级
    df["priority"] = get_class_priority(df["class_code"])

    # 排序并去重，保留每个(qry_gene, ref_gene)对中优先级最高的
    df.sort_values("priority", inplace=True)
    df.drop_duplicates(
        subset=["qry_gene_id", "ref_gene_id"], keep="first", inplace=True
    )

    # 按qry_gene_id和优先级排序，以匹配原始输出顺序
    df.sort_values(["qry_gene_id", "priority"], inplace=True)

    return df[["qry_gene_id", "ref_gene_id", "class_code"]]


def write_mapping_pandas(mapping_df, output_file):
    """使用pandas写入映射结果"""
    mapping_df.to_csv(output_file, sep="\t", index=False, header=True)


def main():
    parser = argparse.ArgumentParser(
        description="从gffcompare的.tmap文件中提取基因ID映射 (v2, pandas优化)"
    )
    parser.add_argument("tmap_file", help="输入的.tmap文件")
    parser.add_argument(
        "-o",
        "--output",
        default="gene_id_map.txt",
        help="输出文件(默认: gene_id_map.txt)",
    )
    parser.add_argument(
        "-c",
        "--classes",
        default="=,c,k,m,j",
        help="有效的class codes，逗号分隔(默认: =,c,k,m,j)",
    )

    args = parser.parse_args()

    # 解析有效的class codes
    valid_classes = set(args.classes.split(","))

    # 解析tmap文件
    print(f"解析tmap文件: {args.tmap_file}", file=sys.stderr)
    mappings_df = parse_tmap_pandas(args.tmap_file, valid_classes)

    # 写入结果
    print(f"找到 {len(mappings_df)} 个基因映射", file=sys.stderr)
    write_mapping_pandas(mappings_df, args.output)
    print(f"结果已写入: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
