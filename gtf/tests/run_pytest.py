#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
使用pytest运行所有测试
"""

import sys
import subprocess


def main():
    """使用pytest运行所有测试"""
    print("=" * 70)
    print("GTF工具测试套件 (pytest版本)")
    print("=" * 70)

    # 运行pytest
    result = subprocess.run(
        [
            sys.executable, "-m", "pytest",
            "test_extract_gene_id_map_pytest.py",
            "test_filter_gene_map_by_chromosome_pytest.py",
            "-v"
        ],
        cwd=__file__.replace('run_pytest.py', '')
    )

    return result.returncode == 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
