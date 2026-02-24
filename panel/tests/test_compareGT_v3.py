"""
compareGT-v3.py 的单元测试。

测试覆盖：
- classify_genotype: 基因型分类逻辑
- compare_samples: 样本比较统计计算
- validate_file: 文件验证
- load_vcf_genotypes: VCF 加载（集成测试）
- main: CLI 端到端测试
"""

from pathlib import Path
from unittest.mock import patch

import polars as pl
import pytest
from click.exceptions import Exit as ClickExit

from panel.compareGT_v3 import (
    classify_genotype,
    compare_samples,
    load_vcf_genotypes,
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
def sample_genotype_df() -> pl.DataFrame:
    """创建测试用基因型 DataFrame。"""
    return pl.DataFrame(
        {
            "CHROM": ["chr1", "chr1", "chr1", "chr1"],
            "POS": [100, 200, 300, 400],
            "REF": ["A", "T", "C", "G"],
            "ALT": ["G", "C", "G", "A"],
            "Sample1": ["A/A", "T/T", "C/G", "A/A"],
            "Sample2": ["A/G", "./.", "C/G", "A/A"],
            "Sample3": ["G/G", "C/C", "C/G", "A/A"],
        }
    )


@pytest.fixture
def compare_list_file(tmp_path: Path) -> Path:
    """创建测试用比较列表文件。"""
    compare_path = tmp_path / "compare.txt"
    compare_path.write_text("Sample1\tSample2\nSample1\tSample3\n")
    return compare_path


# ============== classify_genotype 测试 ==============


class TestClassifyGenotype:
    """测试 classify_genotype 函数。"""

    def test_homozygous_genotype(self) -> None:
        """测试纯合基因型分类。"""
        df = pl.DataFrame({"gt": ["A/A", "G/G", "T/T"]})
        result = df.select(classify_genotype("gt").alias("class"))
        assert result["class"].to_list() == ["HOM", "HOM", "HOM"]

    def test_heterozygous_genotype(self) -> None:
        """测试杂合基因型分类。"""
        df = pl.DataFrame({"gt": ["A/G", "C/T", "G/A"]})
        result = df.select(classify_genotype("gt").alias("class"))
        assert result["class"].to_list() == ["HET", "HET", "HET"]

    def test_missing_genotype(self) -> None:
        """测试缺失基因型分类。"""
        df = pl.DataFrame({"gt": ["./.", ".", None]})
        result = df.select(classify_genotype("gt").alias("class"))
        assert result["class"].to_list() == ["MISS", "MISS", "MISS"]

    def test_phased_genotype(self) -> None:
        """测试 phased 基因型（应在加载时已转换为 unphased）。"""
        # 注意：实际使用中 phased 基因型在 load_vcf_genotypes 中已转换
        # 这里测试转换后的结果
        df = pl.DataFrame({"gt": ["A/A", "A/G"]})
        result = df.select(classify_genotype("gt").alias("class"))
        assert result["class"].to_list() == ["HOM", "HET"]


# ============== compare_samples 测试 ==============


class TestCompareSamples:
    """测试 compare_samples 函数。"""

    def test_basic_comparison(self, sample_genotype_df: pl.DataFrame) -> None:
        """测试基本样本比较。"""
        result = compare_samples(sample_genotype_df, "Sample1", "Sample3")

        assert result[0] == "Sample1"  # 样本 A
        assert result[1] == "Sample3"  # 样本 B
        assert result[2] == 4  # 总位点数
        assert result[3] == 4  # 有效位点（Sample3 无缺失）

    def test_with_missing_data(self, sample_genotype_df: pl.DataFrame) -> None:
        """测试含缺失数据的比较。"""
        result = compare_samples(sample_genotype_df, "Sample1", "Sample2")

        assert result[2] == 4  # 总位点数
        assert result[3] == 3  # 有效位点（Sample2 在 pos 200 缺失）

    def test_identical_samples(self) -> None:
        """测试完全相同的样本。"""
        df = pl.DataFrame(
            {
                "CHROM": ["chr1", "chr1"],
                "POS": [100, 200],
                "REF": ["A", "T"],
                "ALT": ["G", "C"],
                "S1": ["A/A", "T/C"],
                "S2": ["A/A", "T/C"],
            }
        )
        result = compare_samples(df, "S1", "S2")

        assert result[8] == 2  # 整体相似度 = 2
        assert result[9] == 100.0  # 整体相似度% = 100%

    def test_completely_different_samples(self) -> None:
        """测试完全不同的样本。"""
        df = pl.DataFrame(
            {
                "CHROM": ["chr1", "chr1"],
                "POS": [100, 200],
                "REF": ["A", "T"],
                "ALT": ["G", "C"],
                "S1": ["A/A", "T/T"],
                "S2": ["G/G", "C/C"],
            }
        )
        result = compare_samples(df, "S1", "S2")

        assert result[8] == 0  # 整体相似度 = 0
        assert result[14] == 2  # 差异位点数 = 2

    def test_all_missing_returns_zeros(self) -> None:
        """测试全部缺失时返回零值。"""
        df = pl.DataFrame(
            {
                "CHROM": ["chr1"],
                "POS": [100],
                "REF": ["A"],
                "ALT": ["G"],
                "S1": ["./."],
                "S2": ["./."],
            }
        )
        result = compare_samples(df, "S1", "S2")

        assert result[3] == 0  # 有效位点 = 0
        assert result[8] == 0  # 整体相似度 = 0
        assert result[9] == 0.0  # 整体相似度% = 0

    def test_homozygous_similarity(self) -> None:
        """测试纯合位点相似度计算。"""
        df = pl.DataFrame(
            {
                "CHROM": ["chr1", "chr1", "chr1"],
                "POS": [100, 200, 300],
                "REF": ["A", "T", "C"],
                "ALT": ["G", "C", "G"],
                "S1": ["A/A", "T/T", "C/C"],  # 全纯合
                "S2": ["A/A", "C/C", "C/C"],  # 全纯合，pos 200 不同
            }
        )
        result = compare_samples(df, "S1", "S2")

        assert result[4] == 3  # A_纯合 = 3
        assert result[6] == 3  # B_纯合 = 3
        assert result[10] == 2  # 纯合相似度 = 2（pos 100, 300 相同）
        assert round(result[11], 1) == 66.7  # 纯合相似度% ≈ 66.7%

    def test_heterozygous_similarity(self) -> None:
        """测试杂合位点相似度计算。"""
        df = pl.DataFrame(
            {
                "CHROM": ["chr1", "chr1"],
                "POS": [100, 200],
                "REF": ["A", "T"],
                "ALT": ["G", "C"],
                "S1": ["A/G", "T/C"],  # 全杂合
                "S2": ["A/G", "T/T"],  # pos 100 杂合，pos 200 纯合
            }
        )
        result = compare_samples(df, "S1", "S2")

        assert result[5] == 2  # A_杂合 = 2
        assert result[7] == 1  # B_杂合 = 1
        assert result[12] == 1  # 杂合相似度 = 1（pos 100 相同）


# ============== validate_file 测试 ==============


class TestValidateFile:
    """测试 validate_file 函数。"""

    def test_valid_file(self, sample_vcf_file: Path) -> None:
        """测试有效文件通过验证。"""
        # 不应抛出异常
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


# ============== load_vcf_genotypes 测试 ==============


class TestLoadVcfGenotypes:
    """测试 load_vcf_genotypes 函数。"""

    def test_load_basic_vcf(self, sample_vcf_file: Path) -> None:
        """测试加载基本 VCF 文件。"""
        # 使用 patch 禁用 Progress 输出
        with patch("panel.compareGT_v3.Progress"):
            df, samples = load_vcf_genotypes(sample_vcf_file)

        assert samples == ["Sample1", "Sample2", "Sample3"]
        assert df.height == 4  # 4 个变异位点
        assert "CHROM" in df.columns
        assert "POS" in df.columns
        assert "Sample1" in df.columns

    def test_phased_genotypes_converted(self, tmp_path: Path) -> None:
        """测试 phased 基因型被转换为 unphased。"""
        vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t100\t.\tA\tG\t100\tPASS\t.\tGT\t0|1
"""
        vcf_path = tmp_path / "phased.vcf"
        vcf_path.write_text(vcf_content)

        with patch("panel.compareGT_v3.Progress"):
            df, _ = load_vcf_genotypes(vcf_path)

        # 确认 | 被转换为 /
        gt_value = df.filter(pl.col("POS") == 100)["S1"][0]
        assert "|" not in gt_value
        assert "/" in gt_value


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

        # 验证输出内容
        content = output_file.read_text()
        lines = content.strip().split("\n")
        assert len(lines) == 4  # 1 header + 3 pairs (C(3,2) = 3)
        assert "Sample1,Sample2" in content or "Sample2,Sample1" in content

    def test_main_with_compare_list(
        self, sample_vcf_file: Path, compare_list_file: Path, tmp_path: Path
    ) -> None:
        """测试使用比较列表。"""
        from typer.testing import CliRunner

        from panel.compareGT_v3 import app

        runner = CliRunner()
        output_file = tmp_path / "output.csv"

        result = runner.invoke(
            app, [str(sample_vcf_file), str(output_file), "-c", str(compare_list_file)]
        )

        assert result.exit_code == 0

        content = output_file.read_text()
        lines = content.strip().split("\n")
        assert len(lines) == 3  # 1 header + 2 pairs from compare list

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
