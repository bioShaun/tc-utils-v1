"""
compareGT_v3.py 的单元测试。

测试覆盖：
- compute_all_pairs_stats: numpy 向量化统计计算
- format_results: 结果格式化
- validate_file: 文件验证
- load_vcf_as_matrix: VCF 分块加载
- get_vcf_samples: 获取样本列表
- main: CLI 端到端测试
"""

from pathlib import Path

import numpy as np
import pytest
from click.exceptions import Exit as ClickExit

from panel.compareGT_v3 import (
    GT_HET,
    GT_HOM_ALT,
    GT_HOM_REF,
    GT_UNKNOWN,
    compute_all_pairs_stats,
    format_results,
    get_vcf_samples,
    load_vcf_as_matrix,
    validate_file,
)

# ============== Fixtures ==============


@pytest.fixture
def sample_vcf_content() -> str:
    """标准测试 VCF 内容。"""
    return """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3
chr1\t100\t.\tA\tG\t100\tPASS\t.\tGT\t0/0\t0/1\t1/1
chr1\t200\t.\tT\tC\t100\tPASS\t.\tGT\t0/0\t./.\t1/1
chr1\t300\t.\tC\tG\t100\tPASS\t.\tGT\t0/1\t0/1\t0/1
chr1\t400\t.\tG\tA\t100\tPASS\t.\tGT\t1/1\t1/1\t1/1
"""


@pytest.fixture
def sample_vcf_file(tmp_path: Path, sample_vcf_content: str) -> Path:
    """创建测试 VCF 文件。"""
    vcf_path = tmp_path / "test.vcf"
    vcf_path.write_text(sample_vcf_content)
    return vcf_path


@pytest.fixture
def multi_allelic_vcf_content() -> str:
    """包含 multi-allelic 位点的 VCF 内容。"""
    return """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t100\t.\tA\tG\t100\tPASS\t.\tGT\t0/0\t0/1
chr1\t200\t.\tT\tC,G\t100\tPASS\t.\tGT\t0/1\t0/2
chr1\t300\t.\tC\tG\t100\tPASS\t.\tGT\t1/1\t1/1
"""


@pytest.fixture
def multi_allelic_vcf_file(tmp_path: Path, multi_allelic_vcf_content: str) -> Path:
    """创建包含 multi-allelic 位点的测试 VCF 文件。"""
    vcf_path = tmp_path / "multi.vcf"
    vcf_path.write_text(multi_allelic_vcf_content)
    return vcf_path


@pytest.fixture
def compare_list_file(tmp_path: Path) -> Path:
    """创建测试用比较列表文件。"""
    compare_path = tmp_path / "compare.txt"
    compare_path.write_text("Sample1\tSample2\nSample1\tSample3\n")
    return compare_path


# ============== compute_all_pairs_stats 测试 ==============


class TestComputeAllPairsStats:
    """测试 compute_all_pairs_stats 函数。"""

    def test_basic_comparison(self) -> None:
        """测试基本统计计算。"""
        # 创建测试矩阵: 4 variants × 3 samples
        # gt_types: 0=HOM_REF, 1=HET, 2=UNKNOWN, 3=HOM_ALT
        gt_matrix = np.array(
            [
                [GT_HOM_REF, GT_HET, GT_HOM_ALT],  # pos 100: 0/0, 0/1, 1/1
                [GT_HOM_REF, GT_UNKNOWN, GT_HOM_ALT],  # pos 200: 0/0, ./., 1/1
                [GT_HET, GT_HET, GT_HET],  # pos 300: 0/1, 0/1, 0/1
                [GT_HOM_ALT, GT_HOM_ALT, GT_HOM_ALT],  # pos 400: 1/1, 1/1, 1/1
            ],
            dtype=np.int8,
        )

        # 比较 Sample1(idx=0) vs Sample2(idx=1)
        pair_indices = np.array([[0, 1]], dtype=np.int32)
        stats = compute_all_pairs_stats(gt_matrix, pair_indices)

        # Sample1 vs Sample2:
        # pos 100: HOM_REF vs HET -> valid, not equal
        # pos 200: HOM_REF vs UNKNOWN -> invalid (skip)
        # pos 300: HET vs HET -> valid, equal
        # pos 400: HOM_ALT vs HOM_ALT -> valid, equal
        assert stats["non_miss"][0] == 3  # 3 个有效位点
        assert stats["a_hom"][0] == 2  # Sample1 有 2 个纯合
        assert stats["a_het"][0] == 1  # Sample1 有 1 个杂合

    def test_identical_samples(self) -> None:
        """测试完全相同的样本。"""
        gt_matrix = np.array(
            [
                [GT_HOM_REF, GT_HOM_REF],
                [GT_HET, GT_HET],
                [GT_HOM_ALT, GT_HOM_ALT],
            ],
            dtype=np.int8,
        )

        pair_indices = np.array([[0, 1]], dtype=np.int32)
        stats = compute_all_pairs_stats(gt_matrix, pair_indices)

        assert stats["non_miss"][0] == 3
        assert stats["homo_equal"][0] == 2  # 2 个纯合相等
        assert stats["het_equal"][0] == 1  # 1 个杂合相等

    def test_completely_different_samples(self) -> None:
        """测试完全不同的样本。"""
        gt_matrix = np.array(
            [
                [GT_HOM_REF, GT_HOM_ALT],
                [GT_HOM_ALT, GT_HOM_REF],
            ],
            dtype=np.int8,
        )

        pair_indices = np.array([[0, 1]], dtype=np.int32)
        stats = compute_all_pairs_stats(gt_matrix, pair_indices)

        assert stats["non_miss"][0] == 2
        assert stats["homo_equal"][0] == 0  # 无相等的纯合
        assert stats["het_equal"][0] == 0  # 无杂合

    def test_all_missing(self) -> None:
        """测试全部缺失的情况。"""
        gt_matrix = np.array(
            [
                [GT_UNKNOWN, GT_UNKNOWN],
                [GT_UNKNOWN, GT_UNKNOWN],
            ],
            dtype=np.int8,
        )

        pair_indices = np.array([[0, 1]], dtype=np.int32)
        stats = compute_all_pairs_stats(gt_matrix, pair_indices)

        assert stats["non_miss"][0] == 0
        assert stats["homo_equal"][0] == 0
        assert stats["het_equal"][0] == 0

    def test_multiple_pairs(self) -> None:
        """测试多个样本对同时计算。"""
        gt_matrix = np.array(
            [
                [GT_HOM_REF, GT_HET, GT_HOM_ALT],
                [GT_HET, GT_HET, GT_HET],
            ],
            dtype=np.int8,
        )

        # 比较 (0,1), (0,2), (1,2)
        pair_indices = np.array([[0, 1], [0, 2], [1, 2]], dtype=np.int32)
        stats = compute_all_pairs_stats(gt_matrix, pair_indices)

        assert len(stats["non_miss"]) == 3
        assert stats["non_miss"].tolist() == [2, 2, 2]


# ============== format_results 测试 ==============


class TestFormatResults:
    """测试 format_results 函数。"""

    def test_basic_format(self) -> None:
        """测试基本结果格式化。"""
        pairs = [("S1", "S2")]
        # 列: total, non_miss, a_hom, a_het, b_hom, b_het, both_hom, homo_equal, het_sites, het_equal
        accumulators = np.array([[100, 80, 60, 20, 50, 30, 40, 35, 40, 30]], dtype=np.int64)

        results = format_results(pairs, accumulators)

        assert len(results) == 1
        assert results[0][0] == "S1"
        assert results[0][1] == "S2"
        assert results[0][2] == 100  # total_sites
        assert results[0][3] == 80  # non_miss
        assert results[0][8] == 65  # total_equal = homo_equal + het_equal = 35 + 30

    def test_zero_non_miss(self) -> None:
        """测试无有效位点时返回零值。"""
        pairs = [("S1", "S2")]
        accumulators = np.array([[100, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype=np.int64)

        results = format_results(pairs, accumulators)

        assert results[0][3] == 0  # non_miss
        assert results[0][8] == 0  # total_equal
        assert results[0][9] == 0.0  # total_equal_pct

    def test_percentage_calculation(self) -> None:
        """测试百分比计算。"""
        pairs = [("S1", "S2")]
        # 10 个有效位点，8 个相等（5 纯合 + 3 杂合）
        accumulators = np.array([[10, 10, 6, 4, 6, 4, 6, 5, 4, 3]], dtype=np.int64)

        results = format_results(pairs, accumulators)

        assert results[0][9] == 80.0  # 整体相似度% = 8/10 * 100


# ============== validate_file 测试 ==============


class TestValidateFile:
    """测试 validate_file 函数。"""

    def test_valid_file(self, sample_vcf_file: Path) -> None:
        """测试有效文件通过验证。"""
        validate_file(sample_vcf_file, "测试文件")

    def test_nonexistent_file(self, tmp_path: Path) -> None:
        """测试不存在的文件。"""
        fake_path = tmp_path / "nonexistent.vcf"
        with pytest.raises(ClickExit):
            validate_file(fake_path, "测试文件")

    def test_directory_instead_of_file(self, tmp_path: Path) -> None:
        """测试目录而非文件。"""
        with pytest.raises(ClickExit):
            validate_file(tmp_path, "测试文件")


# ============== get_vcf_samples 测试 ==============


class TestGetVcfSamples:
    """测试 get_vcf_samples 函数。"""

    def test_get_samples(self, sample_vcf_file: Path) -> None:
        """测试获取样本列表。"""
        samples = get_vcf_samples(sample_vcf_file)
        assert samples == ["Sample1", "Sample2", "Sample3"]


# ============== load_vcf_as_matrix 测试 ==============


class TestLoadVcfAsMatrix:
    """测试 load_vcf_as_matrix 函数。"""

    def test_load_basic_vcf(self, sample_vcf_file: Path) -> None:
        """测试加载基本 VCF 文件。"""
        chunks = list(load_vcf_as_matrix(sample_vcf_file, chunk_size=100))

        assert len(chunks) == 1  # 只有 4 个变异，一个块
        assert chunks[0].shape == (4, 3)  # 4 variants × 3 samples
        assert chunks[0].dtype == np.int8

    def test_chunk_size(self, sample_vcf_file: Path) -> None:
        """测试分块大小控制。"""
        chunks = list(load_vcf_as_matrix(sample_vcf_file, chunk_size=2))

        assert len(chunks) == 2  # 4 个变异分成 2 块
        assert chunks[0].shape[0] == 2
        assert chunks[1].shape[0] == 2

    def test_gt_types_encoding(self, sample_vcf_file: Path) -> None:
        """测试 gt_types 编码正确性。"""
        chunks = list(load_vcf_as_matrix(sample_vcf_file, chunk_size=100))
        gt_matrix = chunks[0]

        # pos 100: Sample1=0/0, Sample2=0/1, Sample3=1/1
        assert gt_matrix[0, 0] == GT_HOM_REF  # 0/0
        assert gt_matrix[0, 1] == GT_HET  # 0/1
        assert gt_matrix[0, 2] == GT_HOM_ALT  # 1/1

        # pos 200: Sample1=0/0, Sample2=./., Sample3=1/1
        assert gt_matrix[1, 1] == GT_UNKNOWN  # ./.

    def test_multi_allelic_marked_as_unknown(self, multi_allelic_vcf_file: Path) -> None:
        """测试 multi-allelic 位点被标记为 UNKNOWN。"""
        chunks = list(load_vcf_as_matrix(multi_allelic_vcf_file, chunk_size=100))
        gt_matrix = chunks[0]

        # pos 200 是 multi-allelic (ALT=C,G)，应该被标记为 UNKNOWN
        assert gt_matrix[1, 0] == GT_UNKNOWN
        assert gt_matrix[1, 1] == GT_UNKNOWN


# ============== CLI 端到端测试 ==============


class TestCLI:
    """CLI 端到端测试。"""

    def test_main_all_pairs(self, sample_vcf_file: Path, tmp_path: Path) -> None:
        """测试比较所有样本对。"""
        from typer.testing import CliRunner

        from panel.compareGT_v3 import app

        runner = CliRunner()
        output_file = tmp_path / "output.csv"

        result = runner.invoke(app, [str(sample_vcf_file), str(output_file)])

        assert result.exit_code == 0
        assert output_file.exists()

        content = output_file.read_text()
        lines = content.strip().split("\n")
        assert len(lines) == 4  # 1 header + 3 pairs (C(3,2) = 3)
        assert "Sample1,Sample2" in content or "Sample2,Sample1" in content

    def test_main_with_compare_list(self, sample_vcf_file: Path, compare_list_file: Path, tmp_path: Path) -> None:
        """测试使用比较列表。"""
        from typer.testing import CliRunner

        from panel.compareGT_v3 import app

        runner = CliRunner()
        output_file = tmp_path / "output.csv"

        result = runner.invoke(app, [str(sample_vcf_file), str(output_file), "-c", str(compare_list_file)])

        assert result.exit_code == 0

        content = output_file.read_text()
        lines = content.strip().split("\n")
        assert len(lines) == 3  # 1 header + 2 pairs

    def test_main_with_chunk_size(self, sample_vcf_file: Path, tmp_path: Path) -> None:
        """测试自定义分块大小。"""
        from typer.testing import CliRunner

        from panel.compareGT_v3 import app

        runner = CliRunner()
        output_file = tmp_path / "output.csv"

        result = runner.invoke(app, [str(sample_vcf_file), str(output_file), "--chunk-size", "2"])

        assert result.exit_code == 0
        assert output_file.exists()

    def test_main_nonexistent_vcf(self, tmp_path: Path) -> None:
        """测试不存在的 VCF 文件。"""
        from typer.testing import CliRunner

        from panel.compareGT_v3 import app

        runner = CliRunner()
        output_file = tmp_path / "output.csv"
        fake_vcf = tmp_path / "nonexistent.vcf"

        result = runner.invoke(app, [str(fake_vcf), str(output_file)])

        assert result.exit_code == 1

    def test_main_verbose_mode(self, sample_vcf_file: Path, tmp_path: Path) -> None:
        """测试详细日志模式。"""
        from typer.testing import CliRunner

        from panel.compareGT_v3 import app

        runner = CliRunner()
        output_file = tmp_path / "output.csv"

        result = runner.invoke(app, [str(sample_vcf_file), str(output_file), "-v"])

        assert result.exit_code == 0

    def test_multi_allelic_excluded(self, multi_allelic_vcf_file: Path, tmp_path: Path) -> None:
        """测试 multi-allelic 位点被排除出有效位点。"""
        from typer.testing import CliRunner

        from panel.compareGT_v3 import app

        runner = CliRunner()
        output_file = tmp_path / "output.csv"

        result = runner.invoke(app, [str(multi_allelic_vcf_file), str(output_file)])

        assert result.exit_code == 0

        content = output_file.read_text()
        lines = content.strip().split("\n")
        # 解析结果行
        data_line = lines[1].split(",")
        total_sites = int(data_line[2])
        valid_sites = int(data_line[3])

        assert total_sites == 3  # 总位点数保持 3
        assert valid_sites == 2  # 有效位点 = 2（排除 1 个 multi-allelic）
