import gzip
from pathlib import Path

import pytest

from transgene.depth_compare_v2 import (
    DepthResult,
    DepthStats,
    classify_genotype,
    process_single_sample,
    read_depth_file,
)

BAMDST_HEADER = "#Chr\tPos\tRaw Depth"


def create_depth_file(path: Path, depths: list[int]):
    """创建测试用的 bamdst depth.tsv.gz 文件"""
    with gzip.open(path, "wt") as f:
        f.write(BAMDST_HEADER + "\n")
        for i, depth in enumerate(depths, start=1):
            f.write(f"chr1\t{i}\t{depth}\n")


# ── classify_genotype ──────────────────────────────────────────────


class TestClassifyGenotype:
    """测试基因型分类逻辑。"""

    @pytest.mark.parametrize(
        "coverage, ratio, expected",
        [
            (0.0, 0.75, "非转基因"),
            (0.1, 1.4, "非转基因"),
            (0.29, 0.8, "非转基因"),
        ],
    )
    def test_low_coverage_is_non_transgenic(self, coverage, ratio, expected):
        assert classify_genotype(coverage, ratio, tolerance=0.15) == expected

    @pytest.mark.parametrize(
        "ratio, expected",
        [
            (0.61, "转基因杂合"),
            (0.75, "转基因杂合"),
            (0.89, "转基因杂合"),
        ],
    )
    def test_heterozygous(self, ratio, expected):
        assert classify_genotype(0.9, ratio, tolerance=0.15) == expected

    @pytest.mark.parametrize(
        "ratio, expected",
        [
            (1.0, "转基因纯合"),
            (1.2, "转基因纯合"),
            (1.4, "转基因纯合"),
            (1.8, "转基因纯合"),
            (3.0, "转基因纯合"),
        ],
    )
    def test_homozygous(self, ratio, expected):
        assert classify_genotype(0.9, ratio, tolerance=0.15) == expected

    @pytest.mark.parametrize(
        "ratio",
        [0.2, 0.4, 0.95],
    )
    def test_undetermined(self, ratio):
        """杂合上限(0.9)到纯合下限(1.0)之间为未确定"""
        assert classify_genotype(0.9, ratio, tolerance=0.15) == "未确定"

    def test_custom_min_coverage(self):
        """min_coverage=0.5 时，coverage=0.4 应判为非转基因"""
        assert (
            classify_genotype(0.4, 0.75, tolerance=0.15, min_coverage=0.5)
            == "非转基因"
        )

    def test_custom_het_hom_ratio(self):
        """自定义 het_ratio/hom_ratio"""
        assert classify_genotype(0.9, 0.7, tolerance=0.1, het_ratio=0.7) == "转基因杂合"
        assert classify_genotype(0.9, 1.5, tolerance=0.1, hom_ratio=1.5) == "转基因纯合"

    def test_custom_tolerance(self):
        """tolerance=0.05 时边界外的 ratio 应判为未确定"""
        assert classify_genotype(0.9, 0.81, tolerance=0.05) == "未确定"
        assert classify_genotype(0.9, 0.79, tolerance=0.05) == "转基因杂合"

    def test_boundary_coverage_exact(self):
        """coverage 恰好等于 min_coverage 时应进入 ratio 判定"""
        assert classify_genotype(0.3, 0.75, tolerance=0.15) == "转基因杂合"


# ── read_depth_file ────────────────────────────────────────────────


class TestReadDepthFile:
    def test_normal_depths(self, tmp_path):
        """测试正常深度数据的 median 和 coverage"""
        depth_file = tmp_path / "depth.tsv.gz"
        create_depth_file(depth_file, [10, 20, 30, 40, 50])

        result = read_depth_file(depth_file)
        assert result is not None
        assert result.median == 30.0
        assert result.coverage == 1.0

    def test_with_zero_depths(self, tmp_path):
        """测试包含 0 深度的数据"""
        depth_file = tmp_path / "depth.tsv.gz"
        create_depth_file(depth_file, [0, 0, 10, 20, 30, 40, 50, 0])

        result = read_depth_file(depth_file)
        assert result is not None
        assert result.median == 30.0
        assert result.coverage == pytest.approx(5 / 8)

    def test_all_zeros(self, tmp_path):
        """测试全是 0 的数据"""
        depth_file = tmp_path / "depth.tsv.gz"
        create_depth_file(depth_file, [0, 0, 0])

        result = read_depth_file(depth_file)
        assert result is not None
        assert result.median == 0.0
        assert result.coverage == 0.0

    def test_empty_file(self, tmp_path):
        """测试空文件（只有表头）"""
        depth_file = tmp_path / "depth.tsv.gz"
        with gzip.open(depth_file, "wt") as f:
            f.write(BAMDST_HEADER + "\n")

        result = read_depth_file(depth_file)
        assert result is not None
        assert result.median == 0.0
        assert result.coverage == 0.0

    def test_file_not_exists(self, tmp_path):
        result = read_depth_file(tmp_path / "not_exists.tsv.gz")
        assert result is None


# ── process_single_sample ──────────────────────────────────────────


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

        result = process_single_sample(
            transgene_dir, bg_dir, tolerance=0.15, min_coverage=0.3,
            het_ratio=0.75, hom_ratio=1.0,
        )

        assert result is not None
        assert isinstance(result, DepthStats)
        assert result.sample_id == "sample1"
        assert result.transgene_depth == 40.0
        assert result.transgene_coverage == 1.0
        assert result.background_depth == 30.0
        assert abs(result.ratio - 1.333333) < 0.01
        assert result.genotype == "转基因纯合"

    def test_missing_background(self, tmp_path):
        """测试背景文件缺失"""
        transgene_dir = tmp_path / "sample2"
        transgene_dir.mkdir()
        create_depth_file(transgene_dir / "depth.tsv.gz", [10, 20, 30])

        bg_dir = tmp_path / "background"
        bg_dir.mkdir()

        result = process_single_sample(
            transgene_dir, bg_dir, tolerance=0.15, min_coverage=0.3,
            het_ratio=0.75, hom_ratio=1.0,
        )
        assert result is None

    def test_heterozygous_sample(self, tmp_path):
        """测试杂合样本：transgene/background ratio ≈ 0.75"""
        transgene_dir = tmp_path / "sample_het"
        transgene_dir.mkdir()
        create_depth_file(transgene_dir / "depth.tsv.gz", [75, 75, 75, 75, 75])

        bg_dir = tmp_path / "background"
        bg_dir.mkdir()
        bg_sample_dir = bg_dir / "sample_het"
        bg_sample_dir.mkdir()
        create_depth_file(bg_sample_dir / "depth.tsv.gz", [100, 100, 100, 100, 100])

        result = process_single_sample(
            transgene_dir, bg_dir, tolerance=0.15, min_coverage=0.3,
            het_ratio=0.75, hom_ratio=1.0,
        )
        assert result is not None
        assert result.genotype == "转基因杂合"

    def test_homozygous_sample(self, tmp_path):
        """测试纯合样本：transgene/background ratio ≈ 1.4"""
        transgene_dir = tmp_path / "sample_hom"
        transgene_dir.mkdir()
        create_depth_file(transgene_dir / "depth.tsv.gz", [14, 14, 14, 14, 14])

        bg_dir = tmp_path / "background"
        bg_dir.mkdir()
        bg_sample_dir = bg_dir / "sample_hom"
        bg_sample_dir.mkdir()
        create_depth_file(bg_sample_dir / "depth.tsv.gz", [10, 10, 10, 10, 10])

        result = process_single_sample(
            transgene_dir, bg_dir, tolerance=0.15, min_coverage=0.3,
            het_ratio=0.75, hom_ratio=1.0,
        )
        assert result is not None
        assert result.genotype == "转基因纯合"

    def test_non_transgenic_low_coverage(self, tmp_path):
        """测试低覆盖率样本应判定为非转基因"""
        transgene_dir = tmp_path / "sample_low"
        transgene_dir.mkdir()
        # 10 个位点只有 2 个有深度 => coverage = 0.2
        create_depth_file(transgene_dir / "depth.tsv.gz", [0, 0, 0, 0, 0, 0, 0, 0, 10, 10])

        bg_dir = tmp_path / "background"
        bg_dir.mkdir()
        bg_sample_dir = bg_dir / "sample_low"
        bg_sample_dir.mkdir()
        create_depth_file(bg_sample_dir / "depth.tsv.gz", [10, 10, 10, 10, 10, 10, 10, 10, 10, 10])

        result = process_single_sample(
            transgene_dir, bg_dir, tolerance=0.15, min_coverage=0.3,
            het_ratio=0.75, hom_ratio=1.0,
        )
        assert result is not None
        assert result.genotype == "非转基因"
