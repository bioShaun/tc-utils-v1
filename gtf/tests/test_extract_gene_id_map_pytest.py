#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
测试 extract_gene_id_map_from_tmap.py 脚本 (pytest版本)
"""

import os
import sys
import tempfile
import subprocess
import pytest


# 测试数据：tmap文件
TMAP_CONTENT = """ref_gene_id\tclass_code\tqry_gene_id\tqry_gene_name\tqry_chr\tqry_strand\tqry_start\tqry_end\tqry_overlap
TraesCS1A01G000100\t=\tTraesCS1A01G000100\tgene1\tchr1A\t+\t100\t200\t180
TraesCS1A01G000100\tc\tTraesCS1A01G000100LC\tgene1_lc\tchr1A\t+\t150\t250\t150
TraesCS1A01G000200\t=\tTraesCS1A01G000200\tgene2\tchr1B\t+\t300\t400\t180
TraesCS1A01G000200LC\t=\tTraesCS1A01G000200LC\tgene2_lc\tchr1B\t+\t350\t450\t170
TraesCS1A01G000300\tj\tTraesCS1A01G000300\tgene3\tchr2A\t+\t500\t600\t150
TraesCS1A01G000300\tm\tTraesCS1A01G000300\tgene3\tchr2A\t+\t510\t590\t140
TraesCS1A01G000400\t=\tTraesCS1A01G000400\tgene4\tchr3A\t-\t700\t800\t180
TraesCS1A01G000500\t=\tTraesCS1A01G000500\tgene5\tchr4A\t+\t900\t1000\t180
TraesCS1A01G000600\t=\tTraesCS1A01G000600\tgene6\tchr5A\t+\t1100\t1200\t180
"""

# 测试数据：group map文件
# 注意：group_map是ref_gene_id到group_id的映射
GROUP_MAP_CONTENT = """group_id\tgene_id
1\tTraesCS1A01G000100
2\tTraesCS1A01G000200
3\tTraesCS1A01G000200LC
4\tTraesCS1A01G000300
5\tTraesCS1A01G000400
6\tTraesCS1A01G000500
"""

EXPECTED_GENE_MAP = """qry_gene_id\tref_gene_id\tclass_code
TraesCS1A01G000100\tTraesCS1A01G000100\t=
TraesCS1A01G000100LC\tTraesCS1A01G000100\tc
TraesCS1A01G000200\tTraesCS1A01G000200\t=
TraesCS1A01G000200LC\tTraesCS1A01G000200LC\t=
TraesCS1A01G000300\tTraesCS1A01G000300\tm
TraesCS1A01G000400\tTraesCS1A01G000400\t=
TraesCS1A01G000500\tTraesCS1A01G000500\t=
TraesCS1A01G000600\tTraesCS1A01G000600\t=
"""

# 期望的group映射输出：group_id -> qry_gene_id
# 因为group_map建立了ref_gene到group的映射，输出时我们要输出group_id -> qry_gene
EXPECTED_GROUP_MAP = """group_id\tgene_id
1\tTraesCS1A01G000100
1\tTraesCS1A01G000100LC
2\tTraesCS1A01G000200
3\tTraesCS1A01G000200LC
4\tTraesCS1A01G000300
5\tTraesCS1A01G000400
6\tTraesCS1A01G000500
"""


@pytest.fixture
def script_path():
    """获取被测试脚本的路径"""
    return os.path.join(os.path.dirname(__file__), "..", "extract_gene_id_map_from_tmap.py")


@pytest.fixture
def temp_dir():
    """创建临时目录"""
    import tempfile
    import shutil

    tmpdir = tempfile.mkdtemp()
    yield tmpdir
    shutil.rmtree(tmpdir)


def test_without_group_map(temp_dir, script_path):
    """测试不使用group map的情况"""
    tmap_file = os.path.join(temp_dir, "test.tmap")
    output_file = os.path.join(temp_dir, "gene_id_map.txt")

    with open(tmap_file, "w") as f:
        f.write(TMAP_CONTENT)

    # 运行脚本
    result = subprocess.run(
        [sys.executable, script_path, tmap_file, "-o", output_file],
        capture_output=True,
        text=True
    )

    assert result.returncode == 0, f"脚本执行失败: {result.stderr}"
    assert os.path.exists(output_file), "输出文件不存在"

    # 验证内容
    with open(output_file) as f:
        content = f.read()

    lines = [line for line in content.strip().split('\n') if line]
    expected_lines = [line for line in EXPECTED_GENE_MAP.strip().split('\n') if line]

    assert lines == expected_lines, f"输出内容不匹配\nexpected: {expected_lines}\ngot: {lines}"


def test_with_group_map(temp_dir, script_path):
    """测试使用group map的情况"""
    tmap_file = os.path.join(temp_dir, "test.tmap")
    group_map_file = os.path.join(temp_dir, "group_map.txt")
    output_file = os.path.join(temp_dir, "gene_id_map.txt")

    with open(tmap_file, "w") as f:
        f.write(TMAP_CONTENT)

    with open(group_map_file, "w") as f:
        f.write(GROUP_MAP_CONTENT)

    # 运行脚本
    result = subprocess.run(
        [sys.executable, script_path, tmap_file, "-g", group_map_file, "-o", output_file],
        capture_output=True,
        text=True
    )

    assert result.returncode == 0, f"脚本执行失败: {result.stderr}"

    # 验证基因映射输出文件
    gene_output = output_file
    assert os.path.exists(gene_output), "基因映射输出文件不存在"

    # 验证group map输出文件
    group_output = output_file.replace('.txt', '') + '_group_map.txt'
    assert os.path.exists(group_output), f"Group映射输出文件不存在: {group_output}"

    # 验证基因映射内容
    with open(gene_output) as f:
        gene_content = f.read()

    gene_lines = [line for line in gene_content.strip().split('\n') if line]
    expected_gene_lines = [line for line in EXPECTED_GENE_MAP.strip().split('\n') if line]

    assert gene_lines == expected_gene_lines, f"基因映射内容不匹配\nexpected: {expected_gene_lines}\ngot: {gene_lines}"

    # 验证group映射内容
    with open(group_output) as f:
        group_content = f.read()

    group_lines = [line for line in group_content.strip().split('\n') if line]
    expected_group_lines = [line for line in EXPECTED_GROUP_MAP.strip().split('\n') if line]

    assert group_lines == expected_group_lines, f"Group映射内容不匹配\nexpected: {expected_group_lines}\ngot: {group_lines}"


def test_priority_handling(temp_dir, script_path):
    """测试优先级处理（同一对基因多个class code的情况）"""
    # 创建一个有重复基因对的tmap文件
    tmap_with_priority = """ref_gene_id\tclass_code\tqry_gene_id\tqry_gene_name\tqry_chr\tqry_strand\tqry_start\tqry_end\tqry_overlap
TraesCS1A01G000100\t=\tTraesCS1A01G000100\tgene1\tchr1A\t+\t100\t200\t180
TraesCS1A01G000100\tj\tTraesCS1A01G000100\tgene1\tchr1A\t+\t100\t200\t180
TraesCS1A01G000100\tc\tTraesCS1A01G000100LC\tgene1_lc\tchr1A\t+\t150\t250\t150
"""

    tmap_file = os.path.join(temp_dir, "test_priority.tmap")
    output_file = os.path.join(temp_dir, "gene_id_map.txt")

    with open(tmap_file, "w") as f:
        f.write(tmap_with_priority)

    # 运行脚本
    result = subprocess.run(
        [sys.executable, script_path, tmap_file, "-o", output_file],
        capture_output=True,
        text=True
    )

    assert result.returncode == 0, f"脚本执行失败: {result.stderr}"

    # 验证输出文件
    with open(output_file) as f:
        content = f.read()

    # 验证优先级：= (优先级1) 应该保留，j (优先级4) 应该被跳过
    assert "TraesCS1A01G000100\tTraesCS1A01G000100\t=" in content, "高优先级的记录未被保留"
    assert "TraesCS1A01G000100\tTraesCS1A01G000100\tj" not in content, "低优先级的记录被保留了"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
