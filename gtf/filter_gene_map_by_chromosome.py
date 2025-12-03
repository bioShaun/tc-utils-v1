#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
根据染色体信息过滤基因映射

流程:
1. 从GFF文件中提取基因ID到染色体ID的映射
2. 读取genome映射文件（qry_genome, ref_genome）
3. 过滤基因映射，只保留qry_gene和ref_gene染色体都存在于genome映射中的记录
"""

import argparse
import sys
import re


def parse_gff(gff_file):
    """
    解析GFF文件，提取基因ID到染色体ID的映射

    返回: {gene_id: chrom_id}
    """
    gene_to_chrom = {}

    with open(gff_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            chrom = fields[0]
            feature_type = fields[2]

            # 只处理基因或转录本
            if feature_type.lower() in ["gene", "transcript", "mrna"]:
                # 从attributes字段提取gene_id
                attr_str = fields[8]

                # 尝试多种常见的gene_id提取模式
                gene_id = None

                # 模式1: gene_id="xxxxx" 或 gene_id = "xxxxx"
                match = re.search(r'gene_id\s*=\s*"?([^";\s]+)"?', attr_str)
                if match:
                    gene_id = match.group(1)

                # 模式2: gene_id "xxxxx" (没有等号)
                if not gene_id:
                    match = re.search(r'gene_id\s+"?([^";\s]+)"?', attr_str)
                    if match:
                        gene_id = match.group(1)

                # 模式3: ID=xxxxx
                if not gene_id:
                    match = re.search(r'\bID\s*=\s*([^;\s]+)', attr_str)
                    if match:
                        gene_id = match.group(1)

                # 模式4: Name=xxxxx
                if not gene_id:
                    match = re.search(r'\bName\s*=\s*([^;\s]+)', attr_str)
                    if match:
                        gene_id = match.group(1)

                if gene_id:
                    gene_to_chrom[gene_id] = chrom

    return gene_to_chrom


def read_genome_map(genome_map_file):
    """
    读取genome映射文件

    格式: qry_genome,ref_genome (逗号分隔，无表头)
    返回: {qry_genome: ref_genome}
    """
    genome_mapping = {}

    with open(genome_map_file, "r") as f:
        for line in f:
            if not line.strip():
                continue

            fields = line.strip().split(",")
            if len(fields) >= 2:
                qry_genome = fields[0].strip()
                ref_genome = fields[1].strip()
                genome_mapping[qry_genome] = ref_genome

    return genome_mapping


def filter_gene_map(gene_map_file, qry_to_chrom, ref_to_chrom, genome_mapping):
    """
    过滤基因映射，只保留染色体匹配的记录

    返回: [(qry_gene, ref_gene, class_code)]
    """
    filtered = []
    skipped_count = 0
    total_count = 0

    with open(gene_map_file, "r") as f:
        header = f.readline().strip()

        for line in f:
            if not line.strip():
                continue

            total_count += 1
            fields = line.strip().split("\t")

            if len(fields) < 3:
                skipped_count += 1
                continue

            qry_gene = fields[0]
            ref_gene = fields[1]
            class_code = fields[2]

            # 检查基因是否存在于染色体映射中
            if qry_gene not in qry_to_chrom or ref_gene not in ref_to_chrom:
                skipped_count += 1
                continue

            qry_chrom = qry_to_chrom[qry_gene]
            ref_chrom = ref_to_chrom[ref_gene]

            # 检查染色体是否在genome映射中匹配
            if qry_chrom not in genome_mapping:
                skipped_count += 1
                continue

            # 检查query染色体映射到的reference染色体是否匹配当前ref染色体
            if genome_mapping[qry_chrom] != ref_chrom:
                skipped_count += 1
                continue

            # 通过所有检查，保留此记录
            filtered.append((qry_gene, ref_gene, class_code))

    return filtered, total_count, skipped_count


def filter_group_map(group_map_file, gene_map_file, qry_to_chrom, ref_to_chrom, genome_mapping):
    """
    过滤group map，只保留对应的qry_gene通过染色体检查的记录

    注意：group_map文件格式是 group_id -> qry_gene_id（从extract_gene_id_map_from_tmap.py生成）

    返回: [(group_id, qry_gene)]
    """
    # 首先获取通过染色体检查的qry_gene列表
    valid_qry_genes = set()

    with open(gene_map_file, "r") as f:
        header = f.readline().strip()

        for line in f:
            if not line.strip():
                continue

            fields = line.strip().split("\t")
            if len(fields) < 3:
                continue

            qry_gene = fields[0]
            ref_gene = fields[1]

            # 检查基因是否存在于染色体映射中
            if qry_gene not in qry_to_chrom or ref_gene not in ref_to_chrom:
                continue

            qry_chrom = qry_to_chrom[qry_gene]
            ref_chrom = ref_to_chrom[ref_gene]

            # 检查染色体是否在genome映射中匹配
            if qry_chrom not in genome_mapping:
                continue

            # 检查query染色体映射到的reference染色体是否匹配当前ref染色体
            if genome_mapping[qry_chrom] != ref_chrom:
                continue

            # 通过所有检查，添加到有效列表
            valid_qry_genes.add(qry_gene)

    # 过滤group map
    filtered = []
    total_count = 0
    skipped_count = 0

    with open(group_map_file, "r") as f:
        header = f.readline().strip()

        for line in f:
            if not line.strip():
                continue

            total_count += 1
            fields = line.strip().split("\t")

            if len(fields) < 2:
                skipped_count += 1
                continue

            group_id = fields[0]
            qry_gene = fields[1]  # group_map中的gene_id是qry_gene_id

            # 只保留通过染色体检查的qry_gene
            if qry_gene in valid_qry_genes:
                filtered.append((group_id, qry_gene))
            else:
                skipped_count += 1

    return filtered, total_count, skipped_count


def write_filtered_map(filtered_mappings, output_file):
    """写入过滤后的映射结果"""
    with open(output_file, "w") as f:
        f.write("qry_gene_id\tref_gene_id\tclass_code\n")
        for qry_gene, ref_gene, class_code in filtered_mappings:
            f.write(f"{qry_gene}\t{ref_gene}\t{class_code}\n")


def main():
    parser = argparse.ArgumentParser(
        description="根据染色体信息过滤基因映射"
    )
    parser.add_argument(
        "--gene-map",
        required=True,
        help="输入的基因映射文件 (qry_gene_id\\tref_gene_id\\tclass_code)"
    )
    parser.add_argument(
        "--qry-gff",
        required=True,
        help="query基因的GFF文件"
    )
    parser.add_argument(
        "--ref-gff",
        required=True,
        help="reference基因的GFF文件"
    )
    parser.add_argument(
        "--genome-map",
        required=True,
        help="genome映射文件 (qry_genome,ref_genome，无表头)"
    )
    parser.add_argument(
        "--group-map",
        help="输入的group map文件 (group_id\\tgene_id，可选)"
    )
    parser.add_argument(
        "-o",
        "--output",
        help="输出文件(默认: filtered_gene_map.txt)"
    )

    args = parser.parse_args()

    # 设置输出文件
    output_file = args.output if args.output else "filtered_gene_map.txt"

    # 解析文件
    print(f"解析query GFF文件: {args.qry_gff}", file=sys.stderr)
    qry_to_chrom = parse_gff(args.qry_gff)
    print(f"  提取到 {len(qry_to_chrom)} 个基因-染色体映射", file=sys.stderr)

    print(f"解析reference GFF文件: {args.ref_gff}", file=sys.stderr)
    ref_to_chrom = parse_gff(args.ref_gff)
    print(f"  提取到 {len(ref_to_chrom)} 个基因-染色体映射", file=sys.stderr)

    print(f"解析genome映射文件: {args.genome_map}", file=sys.stderr)
    genome_mapping = read_genome_map(args.genome_map)
    print(f"  读取到 {len(genome_mapping)} 个genome映射", file=sys.stderr)

    # 过滤基因映射
    print(f"过滤基因映射文件: {args.gene_map}", file=sys.stderr)
    filtered, total, skipped = filter_gene_map(
        args.gene_map, qry_to_chrom, ref_to_chrom, genome_mapping
    )

    # 写入过滤后的基因映射
    write_filtered_map(filtered, output_file)

    # 如果提供了group map文件，过滤它
    group_output = None
    if args.group_map:
        print(f"过滤group map文件: {args.group_map}", file=sys.stderr)
        group_output = output_file.replace('.txt', '') + '_filtered_group_map.txt'
        filtered_group, group_total, group_skipped = filter_group_map(
            args.group_map, args.gene_map, qry_to_chrom, ref_to_chrom,
            genome_mapping
        )

        # 写入过滤后的group映射
        with open(group_output, "w") as f:
            f.write("group_id\tgene_id\n")
            for group_id, qry_gene in filtered_group:
                f.write(f"{group_id}\t{qry_gene}\n")

    # 输出统计信息
    print(f"\\n=== 基因映射过滤结果 ===", file=sys.stderr)
    print(f"  总记录数: {total}", file=sys.stderr)
    print(f"  跳过记录数: {skipped}", file=sys.stderr)
    print(f"  保留记录数: {len(filtered)}", file=sys.stderr)
    print(f"  过滤率: {skipped/total*100:.2f}%" if total > 0 else "  过滤率: N/A", file=sys.stderr)

    if group_output:
        print(f"\\n=== Group映射过滤结果 ===", file=sys.stderr)
        print(f"  总记录数: {group_total}", file=sys.stderr)
        print(f"  跳过记录数: {group_skipped}", file=sys.stderr)
        print(f"  保留记录数: {len(filtered_group)}", file=sys.stderr)
        print(f"  过滤率: {group_skipped/group_total*100:.2f}%" if group_total > 0 else "  过滤率: N/A", file=sys.stderr)

    print(f"\\n基因映射结果已写入: {output_file}", file=sys.stderr)
    if group_output:
        print(f"Group映射结果已写入: {group_output}", file=sys.stderr)


if __name__ == "__main__":
    main()