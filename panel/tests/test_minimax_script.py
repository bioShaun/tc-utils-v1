#!/usr/bin/env python3
"""
测试 cdsCovEvaluation-bamdst-minimax.py 脚本的功能
"""

import gzip
import sys
import tempfile
from pathlib import Path
import pandas as pd

# 添加路径
sys.path.insert(0, str(Path(__file__).parent / "panel"))

# 直接导入模块文件
import importlib.util
script_path = Path(__file__).parent.parent / "cdsCovEvaluation-bamdst-minimax.py"
spec = importlib.util.spec_from_file_location(
    "cdsCovEvaluation_bamdst_minimax",
    str(script_path)
)
if spec is None or spec.loader is None:
    raise ImportError("Failed to load module spec")
minimax_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(minimax_module)

# 获取类
FileLoader = minimax_module.FileLoader
DataMerger = minimax_module.DataMerger
CoverageAnalyzer = minimax_module.CoverageAnalyzer
CoordinateTransformer = minimax_module.CoordinateTransformer


def create_test_data_region_mode(tmp_dir: Path) -> Path:
    """创建 region.tsv.gz 格式的测试数据"""
    cds_dir = tmp_dir / "cds_cov_region"
    cds_dir.mkdir(parents=True)

    # 创建两个样本目录
    for sample in ["sample1", "sample2"]:
        sample_dir = cds_dir / sample
        sample_dir.mkdir(parents=True)

        # 创建 region.tsv.gz 文件
        region_file = sample_dir / "region.tsv.gz"

        # 写入测试数据：chrom, start, end, coverage
        data = [
            "chr1\t100\t200\t30",  # 样本1：第一区间覆盖度30
            "chr1\t200\t300\t50",  # 样本1：第二区间覆盖度50
            "chr1\t300\t400\t20",  # 样本1：第三区间覆盖度20
            "chr1\t400\t500\t10",  # 样本1：第四区间覆盖度10
        ]

        with gzip.open(region_file, "wt") as f:
            f.write("\n".join(data))

    return cds_dir


def create_test_data_site_mode(tmp_dir: Path) -> Path:
    """创建 depth.tsv.gz 格式的测试数据"""
    cds_dir = tmp_dir / "cds_cov_site"
    cds_dir.mkdir(parents=True)

    # 创建两个样本目录
    for sample in ["sample1", "sample2"]:
        sample_dir = cds_dir / sample
        sample_dir.mkdir(parents=True)

        # 创建 depth.tsv.gz 文件
        depth_file = sample_dir / "depth.tsv.gz"

        # 写入测试数据：#Chr, Pos, Raw Depth, Rmdup depth, Cover depth
        data = [
            "#Chr\tPos\tRaw Depth\tRmdup depth\tCover depth",  # 表头
            "chr1\t101\t35\t32\t30",  # 样本1：位置101，覆盖度30
            "chr1\t201\t55\t52\t50",  # 样本1：位置201，覆盖度50
            "chr1\t301\t25\t22\t20",  # 样本1：位置301，覆盖度20
            "chr1\t401\t15\t12\t10",  # 样本1：位置401，覆盖度10
        ]

        with gzip.open(depth_file, "wt") as f:
            f.write("\n".join(data))

    return cds_dir


def test_region_mode():
    """测试 region.tsv.gz 模式"""
    print("=" * 60)
    print("测试 1: region.tsv.gz 模式 (use_site=False)")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        cds_dir = create_test_data_region_mode(tmp_path)

        # 测试加载
        try:
            bed_df, df_list = FileLoader.load_region_files(cds_dir)
            print("✓ 成功加载 region.tsv.gz 文件")
            print(f"  - BED 区域数: {len(bed_df)}")
            print(f"  - 样本数: {len(df_list)}")
            print(f"  - BED 列: {list(bed_df.columns)}")

            # 合并数据
            merged_df = DataMerger.merge_sample_dataframes(df_list)
            print(f"✓ 成功合并数据")
            print(f"  - 合并后形状: {merged_df.shape}")
            print(f"  - 列名: {list(merged_df.columns)}")

            # 计算统计信息
            stats_df = CoverageAnalyzer.calculate_summary_stats(merged_df)
            print(f"✓ 成功计算统计信息")
            print(f"  - 统计列: {list(stats_df.columns)}")

            # 计算覆盖率
            coverage_ratios = CoverageAnalyzer.calculate_coverage_ratios(
                merged_df, thresholds=[10, 20, 30]
            )
            print(f"✓ 成功计算覆盖率")
            for ratio in coverage_ratios:
                print(f"  - {ratio.name}: {ratio.tolist()}")

            print("\n✅ region.tsv.gz 模式测试通过！\n")

        except Exception as e:
            print(f"\n❌ 测试失败: {e}\n")
            import traceback
            traceback.print_exc()
            sys.exit(1)


def test_site_mode():
    """测试 depth.tsv.gz 模式 (use_site=True)"""
    print("=" * 60)
    print("测试 2: depth.tsv.gz 模式 (use_site=True)")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        cds_dir = create_test_data_site_mode(tmp_path)

        # 测试加载
        try:
            bed_df, df_list = FileLoader.load_depth_files(cds_dir)
            print("✓ 成功加载 depth.tsv.gz 文件")
            print(f"  - BED 区域数: {len(bed_df)}")
            print(f"  - 样本数: {len(df_list)}")
            print(f"  - BED 列: {list(bed_df.columns)}")

            # 打印前几行 BED 数据
            print("\n  BED 数据预览:")
            print(bed_df.to_string(index=False))

            # 合并数据
            merged_df = DataMerger.merge_sample_dataframes(df_list)
            print(f"\n✓ 成功合并数据")
            print(f"  - 合并后形状: {merged_df.shape}")
            print(f"  - 列名: {list(merged_df.columns)}")

            # 打印覆盖矩阵
            print("\n  覆盖矩阵预览:")
            print(merged_df.to_string())

            # 计算统计信息
            stats_df = CoverageAnalyzer.calculate_summary_stats(merged_df)
            print(f"\n✓ 成功计算统计信息")
            print(f"  - 统计列: {list(stats_df.columns)}")
            print("\n  统计信息预览:")
            print(stats_df.to_string())

            # 计算覆盖率
            coverage_ratios = CoverageAnalyzer.calculate_coverage_ratios(
                merged_df, thresholds=[10, 20, 30]
            )
            print(f"\n✓ 成功计算覆盖率")
            for ratio in coverage_ratios:
                print(f"  - {ratio.name}: {ratio.tolist()}")

            print("\n✅ depth.tsv.gz 模式测试通过！\n")

        except Exception as e:
            print(f"\n❌ 测试失败: {e}\n")
            import traceback
            traceback.print_exc()
            sys.exit(1)


def test_coordinate_transformation():
    """测试坐标转换功能"""
    print("=" * 60)
    print("测试 3: 坐标转换功能")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # 创建测试 BED 数据
        bed_df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr2"],
                "start": [100, 200, 300],
                "end": [150, 250, 350],
            }
        )
        print("原始 BED 数据:")
        print(bed_df.to_string(index=False))

        # 创建 split BED 文件
        split_bed_file = tmp_path / "split.bed"
        split_bed_file.write_text(
            "chr1_new\t0\t1000\tchr1\n"
            "chr2_new\t0\t2000\tchr2\n"
        )

        try:
            transformed_df = CoordinateTransformer.merge_chr_coordinates(
                bed_df, split_bed_file
            )
            print("\n✓ 成功转换坐标")
            print("转换后数据:")
            print(transformed_df.to_string(index=False))

            print("\n✅ 坐标转换测试通过！\n")

        except Exception as e:
            print(f"\n❌ 测试失败: {e}\n")
            import traceback
            traceback.print_exc()
            sys.exit(1)


def test_cov_cutoff_filter():
    """测试覆盖率过滤功能"""
    print("=" * 60)
    print("测试 4: 覆盖率过滤功能")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # 创建三个样本，其中一个覆盖率低
        cds_dir = tmp_path / "cds_cov_filter"
        cds_dir.mkdir(parents=True)

        for sample in ["high_cov", "medium_cov", "low_cov"]:
            sample_dir = cds_dir / sample
            sample_dir.mkdir(parents=True)

            # 设置不同覆盖率
            if sample == "low_cov":
                coverage = [5, 8, 3]  # 中位数约6
            elif sample == "medium_cov":
                coverage = [15, 20, 18]  # 中位数约18
            else:  # high_cov
                coverage = [50, 60, 55]  # 中位数约55

            depth_file = sample_dir / "depth.tsv.gz"
            data = [
                "#Chr\tPos\tRaw Depth\tRmdup depth\tCover depth",
                f"chr1\t101\t{coverage[0]+2}\t{coverage[0]+1}\t{coverage[0]}",
                f"chr1\t201\t{coverage[1]+2}\t{coverage[1]+1}\t{coverage[1]}",
                f"chr1\t301\t{coverage[2]+2}\t{coverage[2]+1}\t{coverage[2]}",
            ]

            with gzip.open(depth_file, "wt") as f:
                f.write("\n".join(data))

        try:
            # 测试使用覆盖率阈值 10（应该过滤掉 low_cov）
            bed_df, df_list = FileLoader.load_depth_files(
                cds_dir, cov_cutoff=10
            )
            print(f"✓ 成功应用覆盖率过滤 (阈值: 10)")
            print(f"  - 过滤前样本: 3")
            print(f"  - 过滤后样本: {len(df_list)}")

            if len(df_list) == 2:
                print("  ✓ 成功过滤掉 low_cov 样本")
            else:
                print(f"  ⚠ 预期过滤后样本数为2，实际为{len(df_list)}")

            print("\n✅ 覆盖率过滤测试通过！\n")

        except Exception as e:
            print(f"\n❌ 测试失败: {e}\n")
            import traceback
            traceback.print_exc()
            sys.exit(1)


def main():
    """运行所有测试"""
    print("\n" + "=" * 60)
    print("  cdsCovEvaluation-bamdst-minimax.py 测试套件")
    print("=" * 60 + "\n")

    try:
        # 运行所有测试
        test_region_mode()
        test_site_mode()
        test_coordinate_transformation()
        test_cov_cutoff_filter()

        print("=" * 60)
        print("  🎉 所有测试通过！")
        print("=" * 60)

    except Exception as e:
        print(f"\n💥 测试套件执行失败: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
