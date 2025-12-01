#!/usr/bin/env python3
"""
端到端测试：完整运行 cdsCovEvaluation-bamdst-minimax.py 脚本
"""

import gzip
import tempfile
from pathlib import Path
import subprocess
import sys


def create_test_data(tmp_dir: Path):
    """创建完整的测试数据"""
    cds_dir = tmp_dir / "cds_cov"
    cds_dir.mkdir(parents=True)

    # 创建样本目录和文件
    for sample in ["sample1", "sample2", "sample3"]:
        sample_dir = cds_dir / sample
        sample_dir.mkdir(parents=True)

        # 创建 depth.tsv.gz 文件（use_site=True 模式）
        depth_file = sample_dir / "depth.tsv.gz"
        depth_data = [
            "#Chr\tPos\tRaw Depth\tRmdup depth\tCover depth",
            "chr1\t101\t35\t32\t30",
            "chr1\t201\t55\t52\t50",
            "chr1\t301\t25\t22\t20",
            "chr1\t401\t15\t12\t10",
            "chr2\t101\t45\t42\t40",
            "chr2\t201\t65\t62\t60",
            "chr3\t101\t25\t22\t20",
        ]

        with gzip.open(depth_file, "wt") as f:
            f.write("\n".join(depth_data))

        # 同时创建 region.tsv.gz 文件（use_site=False 模式）
        region_file = sample_dir / "region.tsv.gz"
        region_data = [
            "chr1\t100\t200\t30",
            "chr1\t200\t300\t50",
            "chr1\t300\t400\t20",
            "chr1\t400\t500\t10",
            "chr2\t100\t200\t40",
            "chr2\t200\t300\t60",
            "chr3\t100\t200\t20",
        ]

        with gzip.open(region_file, "wt") as f:
            f.write("\n".join(region_data))

    return cds_dir


def test_end_to_end():
    """端到端测试"""
    print("\n" + "=" * 60)
    print("  端到端测试：完整运行脚本")
    print("=" * 60 + "\n")

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # 创建测试数据
        cds_dir = create_test_data(tmp_path)
        print(f"✓ 创建测试数据: {cds_dir}")

        # 测试 1: use_site=True (depth.tsv.gz 模式)
        print("\n" + "-" * 60)
        print("测试 1: use_site=True (depth.tsv.gz 模式)")
        print("-" * 60)

        out_file1 = tmp_path / "output_site.tsv"

        # 构建脚本路径
        script_path = Path(__file__).parent.parent / "cdsCovEvaluation-bamdst-minimax.py"

        result = subprocess.run(
            [
                sys.executable,
                str(script_path),
                str(cds_dir),
                str(out_file1),
                "--use-site",
                "--cov",
                "10",
            ],
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            print(f"❌ 执行失败 (exit code: {result.returncode})")
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")
            sys.exit(1)

        if not out_file1.exists():
            print("❌ 输出文件未生成")
            sys.exit(1)

        print(f"✓ 成功生成输出文件: {out_file1}")

        # 验证输出文件内容
        import pandas as pd
        df = pd.read_csv(out_file1, sep="\t")
        print(f"✓ 输出文件包含 {len(df)} 行")
        print(f"✓ 列名: {list(df.columns)}")

        # 验证必需列
        required_cols = ["chrom", "start", "end", "min_cov", "max_cov", "mean_cov"]
        for col in required_cols:
            if col not in df.columns:
                print(f"❌ 缺少必需列: {col}")
                sys.exit(1)

        print(f"✓ 包含所有必需列")

        # 测试 2: use_site=False (region.tsv.gz 模式)
        print("\n" + "-" * 60)
        print("测试 2: use_site=False (region.tsv.gz 模式)")
        print("-" * 60)

        out_file2 = tmp_path / "output_region.tsv"

        result = subprocess.run(
            [
                sys.executable,
                str(script_path),
                str(cds_dir),
                str(out_file2),
                "--cov",
                "20",
            ],
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            print(f"❌ 执行失败 (exit code: {result.returncode})")
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")
            sys.exit(1)

        if not out_file2.exists():
            print("❌ 输出文件未生成")
            sys.exit(1)

        print(f"✓ 成功生成输出文件: {out_file2}")

        # 验证输出文件内容
        df2 = pd.read_csv(out_file2, sep="\t")
        print(f"✓ 输出文件包含 {len(df2)} 行")
        print(f"✓ 列名: {list(df2.columns)}")

        # 测试 3: 使用覆盖率过滤
        print("\n" + "-" * 60)
        print("测试 3: 使用覆盖率过滤 (cov_cutoff=15)")
        print("-" * 60)

        out_file3 = tmp_path / "output_filtered.tsv"

        result = subprocess.run(
            [
                sys.executable,
                str(script_path),
                str(cds_dir),
                str(out_file3),
                "--use-site",
                "--cov-cutoff",
                "15",
                "--cov",
                "30",
            ],
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            print(f"❌ 执行失败 (exit code: {result.returncode})")
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")
            sys.exit(1)

        print(f"✓ 成功应用覆盖率过滤")
        print(f"✓ 成功生成输出文件: {out_file3}")

        print("\n" + "=" * 60)
        print("  🎉 端到端测试通过！")
        print("=" * 60)


if __name__ == "__main__":
    test_end_to_end()
