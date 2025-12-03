#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
运行所有GTF工具的测试
"""

import sys
import subprocess


def main():
    """运行所有测试脚本"""
    print("=" * 70)
    print("GTF工具测试套件")
    print("=" * 70)

    test_scripts = [
        "test_extract_gene_id_map.py",
        "test_filter_gene_map_by_chromosome.py",
    ]

    total_passed = 0
    total_failed = 0

    for test_script in test_scripts:
        print(f"\n{'=' * 70}")
        print(f"运行测试: {test_script}")
        print('=' * 70)

        result = subprocess.run(
            [sys.executable, test_script],
            cwd=__file__.replace('run_all_tests.py', '')
        )

        if result.returncode == 0:
            total_passed += 1
            print(f"✓ {test_script} 通过")
        else:
            total_failed += 1
            print(f"✗ {test_script} 失败")

    print("\n" + "=" * 70)
    print("测试总结")
    print("=" * 70)
    print(f"总通过: {total_passed}/{len(test_scripts)}")
    print(f"总失败: {total_failed}/{len(test_scripts)}")
    print("=" * 70)

    return total_failed == 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
