#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
从gffcompare的.tmap文件中提取基因ID映射

规则:
1. class code 为 =, c, k, m, j 认为是同一基因，保留所有配对
2. 可靠度排序: = > c ≈ k > m > j
"""

import argparse
import sys


def get_class_priority(class_code):
    """返回class code的优先级，数值越小优先级越高"""
    priority = {"=": 1, "c": 2, "k": 2, "m": 3, "j": 4}
    return priority.get(class_code, 999)


def parse_tmap(tmap_file, valid_classes={"=", "c", "k", "m", "j"}):
    """
    解析tmap文件，提取基因ID映射

    tmap文件格式(tab分隔):
    ref_gene_id  class_code  qry_gene_id  ...
    """
    gene_mappings = {}

    with open(tmap_file, "r") as f:
        header = f.readline().strip().split("\t")

        # 查找关键列的索引
        try:
            ref_gene_idx = header.index("ref_gene_id")
            class_code_idx = header.index("class_code")
            qry_gene_idx = header.index("qry_gene_id")
        except ValueError:
            # 如果找不到列名，尝试使用标准位置
            # 标准tmap格式: 列1是ref_gene_id, 列3是class_code, 列5是qry_gene_id
            ref_gene_idx = 0
            class_code_idx = 2
            qry_gene_idx = 4

        for line in f:
            if not line.strip():
                continue

            fields = line.strip().split("\t")

            if len(fields) <= max(ref_gene_idx, class_code_idx, qry_gene_idx):
                continue

            ref_gene = fields[ref_gene_idx]
            class_code = fields[class_code_idx]
            qry_gene = fields[qry_gene_idx]

            # 跳过空值和无效class code
            if not ref_gene or ref_gene == "-" or not qry_gene or qry_gene == "-":
                continue

            if class_code not in valid_classes:
                continue

            priority = get_class_priority(class_code)
            pair_key = (qry_gene, ref_gene)

            # 如果配对已存在，保留优先级更高的
            if pair_key not in gene_mappings or priority < gene_mappings[pair_key][1]:
                gene_mappings[pair_key] = (class_code, priority)

    # 转换为列表并按优先级排序
    result = [(qry, ref, cc, pri) for (qry, ref), (cc, pri) in gene_mappings.items()]
    result.sort(key=lambda x: (x[0], x[3]))

    return result


def write_mapping(mappings, output_file):
    """写入映射结果"""
    with open(output_file, "w") as f:
        f.write("qry_gene_id\tref_gene_id\tclass_code\n")
        for qry_gene, ref_gene, class_code, _ in mappings:
            f.write(f"{qry_gene}\t{ref_gene}\t{class_code}\n")


def main():
    parser = argparse.ArgumentParser(
        description="从gffcompare的.tmap文件中提取基因ID映射"
    )
    parser.add_argument("tmap_file", help="输入的.tmap文件")
    parser.add_argument("-o", "--output", help="输出文件(默认: gene_id_map.txt)")
    parser.add_argument(
        "-c",
        "--classes",
        default="=,c,k,m,j",
        help="有效的class codes，逗号分隔(默认: =,c,k,m,j)",
    )

    args = parser.parse_args()

    # 设置输出文件
    output_file = args.output if args.output else "gene_id_map.txt"

    # 解析有效的class codes
    valid_classes = set(args.classes.split(","))

    # 解析tmap文件
    print(f"解析tmap文件: {args.tmap_file}", file=sys.stderr)
    mappings = parse_tmap(args.tmap_file, valid_classes)

    # 写入结果
    print(f"找到 {len(mappings)} 个基因映射", file=sys.stderr)
    write_mapping(mappings, output_file)
    print(f"结果已写入: {output_file}", file=sys.stderr)


if __name__ == "__main__":
    main()
