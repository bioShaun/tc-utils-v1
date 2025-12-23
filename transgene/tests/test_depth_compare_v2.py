import gzip
from pathlib import Path

import polars as pl
import pytest

from transgene.depth_compare_v2 import read_depth_median, process_single_sample


BAMDST_HEADER = "#Chr\tPos\tRaw Depth"


def create_depth_file(path: Path, depths: list[int]):
    """创建测试用的 bamdst depth.tsv.gz 文件"""
    with gzip.open(path, "wt") as f:
        f.write(BAMDST_HEADER + "\n")
        for i, depth in enumerate(depths, start=1):
            f.write(f"chr1\t{i}\t{depth}\n")


class TestReadDepthMedian:
    def test_normal_depths(self, tmp_path):
        """测试正常深度数据"""
        depth_file = tmp_path / "depth.tsv.gz"
        create_depth_file(depth_file, [10, 20, 30, 40, 50])

        result = read_depth_median(depth_file)
        # [10,20,30,40,50], 中位数是 30
        assert result == 30.0

    def test_with_zero_depths(self, tmp_path):
        """测试包含0深度的数据"""
        depth_file = tmp_path / "depth.tsv.gz"
        create_depth_file(depth_file, [0, 0, 10, 20, 30, 40, 50, 0])

        result = read_depth_median(depth_file)
        # 过滤0后 [10,20,30,40,50], 中位数是 30
        assert result == 30.0

    def test_all_zeros(self, tmp_path):
        """测试全是0的数据"""
        depth_file = tmp_path / "depth.tsv.gz"
        create_depth_file(depth_file, [0, 0, 0])

        result = read_depth_median(depth_file)
        assert result == 0.0

    def test_empty_file(self, tmp_path):
        """测试空文件（只有表头）"""
        depth_file = tmp_path / "depth.tsv.gz"
        with gzip.open(depth_file, "wt") as f:
            f.write(BAMDST_HEADER + "\n")

        result = read_depth_median(depth_file)
        assert result == 0.0

    def test_file_not_exists(self, tmp_path):
        """测试文件不存在"""
        result = read_depth_median(tmp_path / "not_exists.tsv.gz")
        assert result is None


class TestProcessSingleSample:
    def test_complete_sample(self, tmp_path):
        """测试完整的样本处理"""
        transgene_dir = tmp_path / "sample1"
        transgene_dir.mkdir()
        create_depth_file(transgene_dir / "depth.tsv.gz", [20, 30, 40, 50, 60])

        bg_dir = tmp_path / "background"
        bg_dir.mkdir()
        bg_sample_dir = bg_dir / "sample1"
        bg_sample_dir.mkdir()
        create_depth_file(bg_sample_dir / "depth.tsv.gz", [10, 20, 30, 40, 50])

        result = process_single_sample(transgene_dir, bg_dir)

        assert result is not None
        assert result["sample_id"] == "sample1"
        # transgene: [20,30,40,50,60] 中位数 40
        assert result["transgene_depth"] == 40.0
        # background: [10,20,30,40,50] 中位数 30
        assert result["background_depth"] == 30.0
        # ratio: 40/30 = 1.333...
        assert abs(result["ratio"] - 1.333333) < 0.01

    def test_missing_background(self, tmp_path):
        """测试背景文件缺失"""
        transgene_dir = tmp_path / "sample2"
        transgene_dir.mkdir()
        create_depth_file(transgene_dir / "depth.tsv.gz", [10, 20, 30])

        bg_dir = tmp_path / "background"
        bg_dir.mkdir()
        # 不创建 sample2 的背景文件

        result = process_single_sample(transgene_dir, bg_dir)
        # 应该返回 None（背景缺失时跳过）
        assert result is None
