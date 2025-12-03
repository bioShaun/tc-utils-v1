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
from pathlib import Path


def get_class_priority(class_code):
    """返回class code的优先级，数值越小优先级越高"""
    priority = {"=": 1, "c": 2, "k": 2, "m": 3, "j": 4}
    return priority.get(class_code, 999)


def parse_tmap(
    tmap_file, valid_classes={"=", "c", "k", "m", "j"}
) -> list[tuple[str, str, str, int]]:
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


def read_group_map(group_map_file):
    """
    读取group map文件，建立ref_gene_id到group_id的映射

    文件格式:
    group_id\tgene_id（ref_gene_id）
    """
    gene_to_group = {}
    with open(group_map_file, "r") as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) >= 2:
                group_id = fields[0]
                gene_id = fields[1]
                gene_to_group[gene_id] = group_id
    return gene_to_group


def write_mapping(
    mappings: list[tuple[str, str, str, int]], output_file: Path, group_map_file=None
):
    """
    写入映射结果
    总是输出qry_gene_id\tref_gene_id\tclass_code格式
    如果提供了group_map_file，额外输出group_id\tgene_id格式
    """
    gene_map_file = f"{output_file.stem}.gene_map.txt"
    # 总是输出原始格式
    with open(gene_map_file, "w", encoding="utf-8") as f:
        f.write("qry_gene_id\tref_gene_id\tclass_code\n")
        for qry_gene, ref_gene, class_code, _ in mappings:
            f.write(f"{qry_gene}\t{ref_gene}\t{class_code}\n")

    # 如果提供了group map文件，额外输出group_id和qry_gene_id的映射
    if group_map_file:
        # 生成group map输出文件名
        group_output = output_file.with_name(output_file.stem + ".group_map.txt")

        gene_to_group = read_group_map(group_map_file)

        with open(group_output, "w", encoding="utf-8") as f:
            f.write("group_id\tgene_id\n")
            for qry_gene, ref_gene, _, _ in mappings:
                # group_map是ref_gene_id到group_id的映射，所以用ref_gene查找
                if ref_gene in gene_to_group:
                    group_id = gene_to_group[ref_gene]
                    f.write(f"{group_id}\t{qry_gene}\n")

        return group_output

    return None


def main():
    parser = argparse.ArgumentParser(
        description="从gffcompare的.tmap文件中提取基因ID映射"
    )
    parser.add_argument("tmap_file", help="输入的.tmap文件")
    parser.add_argument("-o", "--output", help="输出文件(默认: genome_A)")
    parser.add_argument(
        "-g",
        "--group-map",
        help="group map文件，包含group_id到ref_gene_id的映射",
    )
    parser.add_argument(
        "-c",
        "--classes",
        default="=,c,k,m,j",
        help="有效的class codes，逗号分隔(默认: =,c,k,m,j)",
    )

    args = parser.parse_args()

    # 设置输出文件
    output_file = Path(args.output) if args.output else Path("genome_A")
    outdir = output_file.parent
    outdir.mkdir(parents=True, exist_ok=True)

    # 解析有效的class codes
    valid_classes = set(args.classes.split(","))

    # 解析tmap文件
    print(f"解析tmap文件: {args.tmap_file}", file=sys.stderr)
    mappings = parse_tmap(args.tmap_file, valid_classes)

    # 写入结果
    print(f"找到 {len(mappings)} 个基因映射", file=sys.stderr)
    group_output = write_mapping(mappings, output_file, args.group_map)
    print(f"结果已写入: {output_file}", file=sys.stderr)
    if group_output:
        print(f"Group映射已写入: {group_output}", file=sys.stderr)


if __name__ == "__main__":
    main()
