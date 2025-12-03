#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
测试 filter_gene_map_by_chromosome.py 脚本
"""

import os
import sys
import tempfile
import subprocess

# 测试数据：基因映射文件
GENE_MAP_CONTENT = """qry_gene_id\tref_gene_id\tclass_code
Gene1\tRefGene1\t=
Gene2\tRefGene2\t=
Gene3\tRefGene3\tj
Gene4\tRefGene4\t=
Gene5\tRefGene5\t=
Gene6\tRefGene6\t=
Gene7\tRefGene7\t=
Gene8\tRefGene8\t=
"""

# 测试数据：query GFF文件
QRY_GFF_CONTENT = """##gff-version 3
chr1\t.\tgene\t100\t200\t.\t+\t.\tgene_id "Gene1"
chr1\t.\tgene\t300\t400\t.\t+\t.\tgene_id "Gene2"
chr1\t.\tgene\t500\t600\t.\t+\t.\tgene_id "Gene3"
chr2\t.\tgene\t100\t200\t.\t+\t.\tgene_id "Gene4"
chr2\t.\tgene\t300\t400\t.\t+\t.\tgene_id "Gene5"
chr3\t.\tgene\t100\t200\t.\t+\t.\tgene_id "Gene6"
chr3\t.\tgene\t300\t400\t.\t+\t.\tgene_id "Gene7"
chr4\t.\tgene\t100\t200\t.\t+\t.\tgene_id "Gene8"
"""

# 测试数据：ref GFF文件
REF_GFF_CONTENT = """##gff-version 3
RefChr1\t.\tgene\t100\t200\t.\t+\t.\tgene_id "RefGene1"
RefChr1\t.\tgene\t300\t400\t.\t+\t.\tgene_id "RefGene2"
RefChr1\t.\tgene\t500\t600\t.\t+\t.\tgene_id "RefGene3"
RefChr2\t.\tgene\t100\t200\t.\t+\t.\tgene_id "RefGene4"
RefChr2\t.\tgene\t300\t400\t.\t+\t.\tgene_id "RefGene5"
RefChr3\t.\tgene\t100\t200\t.\t+\t.\tgene_id "RefGene6"
RefChr3\t.\tgene\t300\t400\t.\t+\t.\tgene_id "RefGene7"
RefChr4\t.\tgene\t100\t200\t.\t+\t.\tgene_id "RefGene8"
"""

# 测试数据：genome map文件
GENOME_MAP_CONTENT = """chr1,RefChr1
chr2,RefChr2
chr3,RefChr3
"""

# 测试数据：group map文件
# 注意：group_map是ref_gene_id到group_id的映射
GROUP_MAP_CONTENT = """group_id\tgene_id
1\tRefGene1
2\tRefGene2
3\tRefGene3
4\tRefGene4
5\tRefGene5
6\tRefGene6
7\tRefGene7
"""

# 期望的过滤结果：只保留染色体匹配的基因对
EXPECTED_FILTERED_GENE_MAP = """qry_gene_id\tref_gene_id\tclass_code
Gene1\tRefGene1\t=
Gene2\tRefGene2\t=
Gene3\tRefGene3\tj
Gene4\tRefGene4\t=
Gene5\tRefGene5\t=
Gene6\tRefGene6\t=
Gene7\tRefGene7\t=
"""

# 期望的group映射输出：group_id -> qry_gene_id
EXPECTED_FILTERED_GROUP_MAP = """group_id\tgene_id
1\tGene1
2\tGene2
3\tGene3
4\tGene4
5\tGene5
6\tGene6
7\tGene7
"""


def test_filter_without_group_map():
    """测试不使用group map的过滤"""
    print("测试1：不使用group map的过滤...")

    with tempfile.TemporaryDirectory() as tmpdir:
        gene_map_file = os.path.join(tmpdir, "gene_map.txt")
        qry_gff_file = os.path.join(tmpdir, "qry.gff")
        ref_gff_file = os.path.join(tmpdir, "ref.gff")
        genome_map_file = os.path.join(tmpdir, "genome_map.txt")
        output_file = os.path.join(tmpdir, "filtered_gene_map.txt")

        with open(gene_map_file, "w") as f:
            f.write(GENE_MAP_CONTENT)

        with open(qry_gff_file, "w") as f:
            f.write(QRY_GFF_CONTENT)

        with open(ref_gff_file, "w") as f:
            f.write(REF_GFF_CONTENT)

        with open(genome_map_file, "w") as f:
            f.write(GENOME_MAP_CONTENT)

        # 运行脚本
        script_path = os.path.join(os.path.dirname(__file__), "..", "filter_gene_map_by_chromosome.py")
        result = subprocess.run(
            [
                sys.executable, script_path,
                "--gene-map", gene_map_file,
                "--qry-gff", qry_gff_file,
                "--ref-gff", ref_gff_file,
                "--genome-map", genome_map_file,
                "-o", output_file
            ],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            print(f"  ❌ 脚本执行失败")
            print(f"  错误输出: {result.stderr}")
            return False

        # 验证输出文件
        if not os.path.exists(output_file):
            print(f"  ❌ 输出文件不存在")
            return False

        with open(output_file) as f:
            content = f.read()

        # 验证内容
        lines = [line for line in content.strip().split('\n') if line]
        expected_lines = [line for line in EXPECTED_FILTERED_GENE_MAP.strip().split('\n') if line]

        if lines != expected_lines:
            print(f"  ❌ 输出内容不匹配")
            print(f"  期望: {expected_lines}")
            print(f"  实际: {lines}")
            return False

        # 验证统计信息
        if "保留记录数: 7" not in result.stderr:
            print(f"  ❌ 统计信息不正确")
            print(f"  实际输出: {result.stderr}")
            return False

        print("  ✓ 通过")
        return True


def test_filter_with_group_map():
    """测试使用group map的过滤"""
    print("测试2：使用group map的过滤...")

    with tempfile.TemporaryDirectory() as tmpdir:
        gene_map_file = os.path.join(tmpdir, "gene_map.txt")
        qry_gff_file = os.path.join(tmpdir, "qry.gff")
        ref_gff_file = os.path.join(tmpdir, "ref.gff")
        genome_map_file = os.path.join(tmpdir, "genome_map.txt")
        group_map_file = os.path.join(tmpdir, "group_map.txt")
        output_file = os.path.join(tmpdir, "filtered_gene_map.txt")

        with open(gene_map_file, "w") as f:
            f.write(GENE_MAP_CONTENT)

        with open(qry_gff_file, "w") as f:
            f.write(QRY_GFF_CONTENT)

        with open(ref_gff_file, "w") as f:
            f.write(REF_GFF_CONTENT)

        with open(genome_map_file, "w") as f:
            f.write(GENOME_MAP_CONTENT)

        with open(group_map_file, "w") as f:
            f.write(GROUP_MAP_CONTENT)

        # 运行脚本
        script_path = os.path.join(os.path.dirname(__file__), "..", "filter_gene_map_by_chromosome.py")
        result = subprocess.run(
            [
                sys.executable, script_path,
                "--gene-map", gene_map_file,
                "--qry-gff", qry_gff_file,
                "--ref-gff", ref_gff_file,
                "--genome-map", genome_map_file,
                "--group-map", group_map_file,
                "-o", output_file
            ],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            print(f"  ❌ 脚本执行失败")
            print(f"  错误输出: {result.stderr}")
            return False

        # 验证基因映射输出文件
        if not os.path.exists(output_file):
            print(f"  ❌ 基因映射输出文件不存在")
            return False

        # 验证group映射输出文件
        group_output = output_file.replace('.txt', '') + '_filtered_group_map.txt'
        if not os.path.exists(group_output):
            print(f"  ❌ Group映射输出文件不存在: {group_output}")
            return False

        # 验证基因映射内容
        with open(output_file) as f:
            gene_content = f.read()

        gene_lines = [line for line in gene_content.strip().split('\n') if line]
        expected_gene_lines = [line for line in EXPECTED_FILTERED_GENE_MAP.strip().split('\n') if line]

        if gene_lines != expected_gene_lines:
            print(f"  ❌ 基因映射内容不匹配")
            print(f"  期望: {expected_gene_lines}")
            print(f"  实际: {gene_lines}")
            return False

        # 验证group映射内容
        with open(group_output) as f:
            group_content = f.read()

        group_lines = [line for line in group_content.strip().split('\n') if line]
        expected_group_lines = [line for line in EXPECTED_FILTERED_GROUP_MAP.strip().split('\n') if line]

        if group_lines != expected_group_lines:
            print(f"  ❌ Group映射内容不匹配")
            print(f"  期望: {expected_group_lines}")
            print(f"  实际: {group_lines}")
            return False

        # 验证统计信息
        if "基因映射过滤结果" not in result.stderr or "Group映射过滤结果" not in result.stderr:
            print(f"  ❌ 缺少统计信息")
            print(f"  实际输出: {result.stderr}")
            return False

        print("  ✓ 通过")
        return True


def test_chromosome_mismatch():
    """测试染色体不匹配的情况（Gene8的chr4不在genome map中）"""
    print("测试3：染色体不匹配的情况...")

    with tempfile.TemporaryDirectory() as tmpdir:
        gene_map_file = os.path.join(tmpdir, "gene_map.txt")
        qry_gff_file = os.path.join(tmpdir, "qry.gff")
        ref_gff_file = os.path.join(tmpdir, "ref.gff")
        genome_map_file = os.path.join(tmpdir, "genome_map.txt")
        output_file = os.path.join(tmpdir, "filtered_gene_map.txt")

        with open(gene_map_file, "w") as f:
            f.write(GENE_MAP_CONTENT)

        with open(qry_gff_file, "w") as f:
            f.write(QRY_GFF_CONTENT)

        with open(ref_gff_file, "w") as f:
            f.write(REF_GFF_CONTENT)

        with open(genome_map_file, "w") as f:
            f.write(GENOME_MAP_CONTENT)

        # 运行脚本
        script_path = os.path.join(os.path.dirname(__file__), "..", "filter_gene_map_by_chromosome.py")
        result = subprocess.run(
            [
                sys.executable, script_path,
                "--gene-map", gene_map_file,
                "--qry-gff", qry_gff_file,
                "--ref-gff", ref_gff_file,
                "--genome-map", genome_map_file,
                "-o", output_file
            ],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            print(f"  ❌ 脚本执行失败")
            print(f"  错误输出: {result.stderr}")
            return False

        # 验证Gene8被过滤掉
        with open(output_file) as f:
            content = f.read()

        if "Gene8" in content:
            print(f"  ❌ Gene8应该被过滤掉（染色体chr4不在genome map中）")
            return False

        # 验证保留了7条记录（过滤掉Gene8）
        if "保留记录数: 7" not in result.stderr:
            print(f"  ❌ 记录数不正确，应该保留7条")
            print(f"  实际输出: {result.stderr}")
            return False

        print("  ✓ 通过")
        return True


def test_missing_gene_in_gff():
    """测试基因在GFF中不存在的情况"""
    print("测试4：基因在GFF中不存在...")

    gene_map_missing = """qry_gene_id\tref_gene_id\tclass_code
Gene1\tRefGene1\t=
Gene999\tRefGene999\t=
"""

    with tempfile.TemporaryDirectory() as tmpdir:
        gene_map_file = os.path.join(tmpdir, "gene_map.txt")
        qry_gff_file = os.path.join(tmpdir, "qry.gff")
        ref_gff_file = os.path.join(tmpdir, "ref.gff")
        genome_map_file = os.path.join(tmpdir, "genome_map.txt")
        output_file = os.path.join(tmpdir, "filtered_gene_map.txt")

        with open(gene_map_file, "w") as f:
            f.write(gene_map_missing)

        with open(qry_gff_file, "w") as f:
            f.write(QRY_GFF_CONTENT)

        with open(ref_gff_file, "w") as f:
            f.write(REF_GFF_CONTENT)

        with open(genome_map_file, "w") as f:
            f.write(GENOME_MAP_CONTENT)

        # 运行脚本
        script_path = os.path.join(os.path.dirname(__file__), "..", "filter_gene_map_by_chromosome.py")
        result = subprocess.run(
            [
                sys.executable, script_path,
                "--gene-map", gene_map_file,
                "--qry-gff", qry_gff_file,
                "--ref-gff", ref_gff_file,
                "--genome-map", genome_map_file,
                "-o", output_file
            ],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            print(f"  ❌ 脚本执行失败")
            print(f"  错误输出: {result.stderr}")
            return False

        # 验证只保留存在的基因
        with open(output_file) as f:
            content = f.read()

        if "Gene999" in content or "RefGene999" in content:
            print(f"  ❌ 不存在的基因Gene999不应该被保留")
            return False

        if "Gene1\tRefGene1\t=" not in content:
            print(f"  ❌ 存在的基因Gene1应该被保留")
            return False

        print("  ✓ 通过")
        return True


def main():
    print("=" * 60)
    print("测试 filter_gene_map_by_chromosome.py")
    print("=" * 60)

    tests = [
        test_filter_without_group_map,
        test_filter_with_group_map,
        test_chromosome_mismatch,
        test_missing_gene_in_gff,
    ]

    passed = 0
    failed = 0

    for test in tests:
        try:
            if test():
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"  ❌ 测试异常: {e}")
            import traceback
            traceback.print_exc()
            failed += 1

    print("=" * 60)
    print(f"测试结果: {passed} 通过, {failed} 失败")
    print("=" * 60)

    return failed == 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
