#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
测试 filter_gene_map_by_chromosome.py 脚本 (pytest版本)
"""

import os
import sys
import tempfile
import subprocess
import pytest


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
GROUP_MAP_CONTENT = """group_id\tgene_id
1\tGene1
2\tGene2
3\tGene3
4\tGene4
5\tGene5
6\tGene6
7\tGene7
8\tGene8
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

EXPECTED_FILTERED_GROUP_MAP = """group_id\tgene_id
1\tGene1
2\tGene2
3\tGene3
4\tGene4
5\tGene5
6\tGene6
7\tGene7
"""


@pytest.fixture
def script_path():
    """获取被测试脚本的路径"""
    return os.path.join(os.path.dirname(__file__), "..", "filter_gene_map_by_chromosome.py")


@pytest.fixture
def temp_dir():
    """创建临时目录"""
    import tempfile
    import shutil

    tmpdir = tempfile.mkdtemp()
    yield tmpdir
    shutil.rmtree(tmpdir)


def test_filter_without_group_map(temp_dir, script_path):
    """测试不使用group map的过滤"""
    gene_map_file = os.path.join(temp_dir, "gene_map.txt")
    qry_gff_file = os.path.join(temp_dir, "qry.gff")
    ref_gff_file = os.path.join(temp_dir, "ref.gff")
    genome_map_file = os.path.join(temp_dir, "genome_map.txt")
    output_file = os.path.join(temp_dir, "filtered_gene_map.txt")

    with open(gene_map_file, "w") as f:
        f.write(GENE_MAP_CONTENT)

    with open(qry_gff_file, "w") as f:
        f.write(QRY_GFF_CONTENT)

    with open(ref_gff_file, "w") as f:
        f.write(REF_GFF_CONTENT)

    with open(genome_map_file, "w") as f:
        f.write(GENOME_MAP_CONTENT)

    # 运行脚本
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

    assert result.returncode == 0, f"脚本执行失败: {result.stderr}"
    assert os.path.exists(output_file), "输出文件不存在"

    # 验证内容
    with open(output_file) as f:
        content = f.read()

    lines = [line for line in content.strip().split('\n') if line]
    expected_lines = [line for line in EXPECTED_FILTERED_GENE_MAP.strip().split('\n') if line]

    assert lines == expected_lines, f"输出内容不匹配\nexpected: {expected_lines}\ngot: {lines}"

    # 验证统计信息
    assert "保留记录数: 7" in result.stderr, f"统计信息不正确: {result.stderr}"


def test_filter_with_group_map(temp_dir, script_path):
    """测试使用group map的过滤"""
    gene_map_file = os.path.join(temp_dir, "gene_map.txt")
    qry_gff_file = os.path.join(temp_dir, "qry.gff")
    ref_gff_file = os.path.join(temp_dir, "ref.gff")
    genome_map_file = os.path.join(temp_dir, "genome_map.txt")
    group_map_file = os.path.join(temp_dir, "group_map.txt")
    output_file = os.path.join(temp_dir, "filtered_gene_map.txt")

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

    assert result.returncode == 0, f"脚本执行失败: {result.stderr}"

    # 验证基因映射输出文件
    assert os.path.exists(output_file), "基因映射输出文件不存在"

    # 验证group映射输出文件
    group_output = output_file.replace('.txt', '') + '_filtered_group_map.txt'
    assert os.path.exists(group_output), f"Group映射输出文件不存在: {group_output}"

    # 验证基因映射内容
    with open(output_file) as f:
        gene_content = f.read()

    gene_lines = [line for line in gene_content.strip().split('\n') if line]
    expected_gene_lines = [line for line in EXPECTED_FILTERED_GENE_MAP.strip().split('\n') if line]

    assert gene_lines == expected_gene_lines, f"基因映射内容不匹配\nexpected: {expected_gene_lines}\ngot: {gene_lines}"

    # 验证group映射内容
    with open(group_output) as f:
        group_content = f.read()

    group_lines = [line for line in group_content.strip().split('\n') if line]
    expected_group_lines = [line for line in EXPECTED_FILTERED_GROUP_MAP.strip().split('\n') if line]

    assert group_lines == expected_group_lines, f"Group映射内容不匹配\nexpected: {expected_group_lines}\ngot: {group_lines}"

    # 验证统计信息
    assert "基因映射过滤结果" in result.stderr, f"缺少基因映射统计信息: {result.stderr}"
    assert "Group映射过滤结果" in result.stderr, f"缺少Group映射统计信息: {result.stderr}"


def test_chromosome_mismatch(temp_dir, script_path):
    """测试染色体不匹配的情况（Gene8的chr4不在genome map中）"""
    gene_map_file = os.path.join(temp_dir, "gene_map.txt")
    qry_gff_file = os.path.join(temp_dir, "qry.gff")
    ref_gff_file = os.path.join(temp_dir, "ref.gff")
    genome_map_file = os.path.join(temp_dir, "genome_map.txt")
    output_file = os.path.join(temp_dir, "filtered_gene_map.txt")

    with open(gene_map_file, "w") as f:
        f.write(GENE_MAP_CONTENT)

    with open(qry_gff_file, "w") as f:
        f.write(QRY_GFF_CONTENT)

    with open(ref_gff_file, "w") as f:
        f.write(REF_GFF_CONTENT)

    with open(genome_map_file, "w") as f:
        f.write(GENOME_MAP_CONTENT)

    # 运行脚本
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

    assert result.returncode == 0, f"脚本执行失败: {result.stderr}"

    # 验证Gene8被过滤掉
    with open(output_file) as f:
        content = f.read()

    assert "Gene8" not in content, "Gene8应该被过滤掉（染色体chr4不在genome map中）"

    # 验证保留了7条记录（过滤掉Gene8）
    assert "保留记录数: 7" in result.stderr, f"记录数不正确，应该保留7条: {result.stderr}"


def test_missing_gene_in_gff(temp_dir, script_path):
    """测试基因在GFF中不存在的情况"""
    gene_map_missing = """qry_gene_id\tref_gene_id\tclass_code
Gene1\tRefGene1\t=
Gene999\tRefGene999\t=
"""

    gene_map_file = os.path.join(temp_dir, "gene_map.txt")
    qry_gff_file = os.path.join(temp_dir, "qry.gff")
    ref_gff_file = os.path.join(temp_dir, "ref.gff")
    genome_map_file = os.path.join(temp_dir, "genome_map.txt")
    output_file = os.path.join(temp_dir, "filtered_gene_map.txt")

    with open(gene_map_file, "w") as f:
        f.write(gene_map_missing)

    with open(qry_gff_file, "w") as f:
        f.write(QRY_GFF_CONTENT)

    with open(ref_gff_file, "w") as f:
        f.write(REF_GFF_CONTENT)

    with open(genome_map_file, "w") as f:
        f.write(GENOME_MAP_CONTENT)

    # 运行脚本
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

    assert result.returncode == 0, f"脚本执行失败: {result.stderr}"

    # 验证只保留存在的基因
    with open(output_file) as f:
        content = f.read()

    assert "Gene999" not in content, "不存在的基因Gene999不应该被保留"
    assert "Gene1\tRefGene1\t=" in content, "存在的基因Gene1应该被保留"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
