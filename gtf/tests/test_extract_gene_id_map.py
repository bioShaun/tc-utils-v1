#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
测试 extract_gene_id_map_from_tmap.py 脚本
"""

import os
import sys
import tempfile
import subprocess

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
GROUP_MAP_CONTENT = """group_id\tgene_id
1\tTraesCS1A01G000100
2\tTraesCS1A01G000100LC
3\tTraesCS1A01G000200
4\tTraesCS1A01G000200LC
5\tTraesCS1A01G000300
6\tTraesCS1A01G000400
7\tTraesCS1A01G000500
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

EXPECTED_GROUP_MAP = """group_id\tgene_id
1\tTraesCS1A01G000100
2\tTraesCS1A01G000100LC
3\tTraesCS1A01G000200
4\tTraesCS1A01G000200LC
5\tTraesCS1A01G000300
6\tTraesCS1A01G000400
7\tTraesCS1A01G000500
"""


def test_without_group_map():
    """测试不使用group map的情况"""
    print("测试1：不使用group map...")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmap_file = os.path.join(tmpdir, "test.tmap")
        output_file = os.path.join(tmpdir, "gene_id_map.txt")

        with open(tmap_file, "w") as f:
            f.write(TMAP_CONTENT)

        # 运行脚本
        script_path = os.path.join(os.path.dirname(__file__), "..", "extract_gene_id_map_from_tmap.py")
        result = subprocess.run(
            [sys.executable, script_path, tmap_file, "-o", output_file],
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

        # 验证内容（去除空行和顺序无关的差异）
        lines = [line for line in content.strip().split('\n') if line]
        expected_lines = [line for line in EXPECTED_GENE_MAP.strip().split('\n') if line]

        if lines != expected_lines:
            print(f"  ❌ 输出内容不匹配")
            print(f"  期望: {expected_lines}")
            print(f"  实际: {lines}")
            return False

        print("  ✓ 通过")
        return True


def test_with_group_map():
    """测试使用group map的情况"""
    print("测试2：使用group map...")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmap_file = os.path.join(tmpdir, "test.tmap")
        group_map_file = os.path.join(tmpdir, "group_map.txt")
        output_file = os.path.join(tmpdir, "gene_id_map.txt")

        with open(tmap_file, "w") as f:
            f.write(TMAP_CONTENT)

        with open(group_map_file, "w") as f:
            f.write(GROUP_MAP_CONTENT)

        # 运行脚本
        script_path = os.path.join(os.path.dirname(__file__), "..", "extract_gene_id_map_from_tmap.py")
        result = subprocess.run(
            [sys.executable, script_path, tmap_file, "-g", group_map_file, "-o", output_file],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            print(f"  ❌ 脚本执行失败")
            print(f"  错误输出: {result.stderr}")
            return False

        # 验证主输出文件
        gene_output = output_file
        if not os.path.exists(gene_output):
            print(f"  ❌ 基因映射输出文件不存在")
            return False

        # 验证group map输出文件
        group_output = output_file.replace('.txt', '') + '_group_map.txt'
        if not os.path.exists(group_output):
            print(f"  ❌ Group映射输出文件不存在: {group_output}")
            return False

        # 验证基因映射内容
        with open(gene_output) as f:
            gene_content = f.read()

        gene_lines = [line for line in gene_content.strip().split('\n') if line]
        expected_gene_lines = [line for line in EXPECTED_GENE_MAP.strip().split('\n') if line]

        if gene_lines != expected_gene_lines:
            print(f"  ❌ 基因映射内容不匹配")
            print(f"  期望: {expected_gene_lines}")
            print(f"  实际: {gene_lines}")
            return False

        # 验证group映射内容
        with open(group_output) as f:
            group_content = f.read()

        group_lines = [line for line in group_content.strip().split('\n') if line]
        expected_group_lines = [line for line in EXPECTED_GROUP_MAP.strip().split('\n') if line]

        if group_lines != expected_group_lines:
            print(f"  ❌ Group映射内容不匹配")
            print(f"  期望: {expected_group_lines}")
            print(f"  实际: {group_lines}")
            return False

        print("  ✓ 通过")
        return True


def test_priority_handling():
    """测试优先级处理（同一对基因多个class code的情况）"""
    print("测试3：优先级处理...")

    # 创建一个有重复基因对的tmap文件
    tmap_with_priority = """ref_gene_id\tclass_code\tqry_gene_id\tqry_gene_name\tqry_chr\tqry_strand\tqry_start\tqry_end\tqry_overlap
TraesCS1A01G000100\t=\tTraesCS1A01G000100\tgene1\tchr1A\t+\t100\t200\t180
TraesCS1A01G000100\tj\tTraesCS1A01G000100\tgene1\tchr1A\t+\t100\t200\t180
TraesCS1A01G000100\tc\tTraesCS1A01G000100LC\tgene1_lc\tchr1A\t+\t150\t250\t150
"""

    with tempfile.TemporaryDirectory() as tmpdir:
        tmap_file = os.path.join(tmpdir, "test_priority.tmap")
        output_file = os.path.join(tmpdir, "gene_id_map.txt")

        with open(tmap_file, "w") as f:
            f.write(tmap_with_priority)

        # 运行脚本
        script_path = os.path.join(os.path.dirname(__file__), "..", "extract_gene_id_map_from_tmap.py")
        result = subprocess.run(
            [sys.executable, script_path, tmap_file, "-o", output_file],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            print(f"  ❌ 脚本执行失败")
            print(f"  错误输出: {result.stderr}")
            return False

        # 验证输出文件
        with open(output_file) as f:
            content = f.read()

        # 验证优先级：= (优先级1) 应该保留，j (优先级4) 应该被跳过
        if "TraesCS1A01G000100\tTraesCS1A01G000100\t=" in content:
            if "TraesCS1A01G000100\tTraesCS1A01G000100\tj" in content:
                print(f"  ❌ 优先级处理不正确，低优先级的记录被保留了")
                return False
        else:
            print(f"  ❌ 高优先级的记录未被保留")
            return False

        print("  ✓ 通过")
        return True


def main():
    print("=" * 60)
    print("测试 extract_gene_id_map_from_tmap.py")
    print("=" * 60)

    tests = [
        test_without_group_map,
        test_with_group_map,
        test_priority_handling,
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
