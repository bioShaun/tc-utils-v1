#!/usr/bin/env python3
"""
VCF分析工具的测试模块

此模块包含对vcf_analyzer.py中函数的单元测试。
"""

import csv
import io
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

# 导入被测试的模块
import panel.simple_vcf_stats_large as va


class MockVariant:
    """模拟cyvcf2.Variant对象"""

    def __init__(self, chrom, pos, ref, alt, gt_types, aaf):
        self.CHROM = chrom
        self.POS = pos
        self.REF = ref
        self.ALT = alt
        self.gt_types = np.array(gt_types)
        self.aaf = aaf
        self.is_snp = True


@pytest.fixture
def mock_variant_snp():
    """创建模拟SNP变异位点"""
    return MockVariant("chr1", 1000, "A", ["T"], [0, 1, 1, 2, 3], 0.4)


@pytest.fixture
def mock_variant_non_snp():
    """创建模拟非SNP变异位点"""
    variant = MockVariant("chr1", 2000, "A", ["TG"], [0, 1, 2, 0, 1], 0.3)
    variant.is_snp = False
    return variant


@pytest.fixture
def mock_variant_multi_alt():
    """创建模拟多等位基因变异位点"""
    return MockVariant("chr1", 3000, "A", ["T", "G"], [0, 1, 2, 0, 1], 0.3)


@pytest.fixture
def mock_variants():
    """创建模拟变异位点列表"""
    return [
        # 正常SNP
        MockVariant("chr1", 1000, "A", ["T"], [0, 1, 1, 2, 3], 0.4),
        # 非SNP (将被跳过)
        MockVariant("chr1", 2000, "A", ["TG"], [0, 1, 2, 0, 1], 0.3),
        # 多个替代等位基因 (将被跳过)
        MockVariant("chr1", 3000, "A", ["T", "G"], [0, 1, 2, 0, 1], 0.3),
        # 正常SNP，有缺失值
        MockVariant("chr2", 1500, "C", ["G"], [0, 3, 1, 2, 3], 0.25),
    ]


def test_vcf_stats_to_dict():
    """测试VcfStats转换为字典"""
    stats = va.VcfStats(
        chrom="chr1", pos=1000, alleles="A,T", missing=0.05, het=0.3, maf=0.25
    )

    data = stats.to_dict()
    assert data["chrom"] == "chr1"
    assert data["pos"] == 1000
    assert data["alleles"] == "A,T"
    assert data["missing"] == 0.05
    assert data["het"] == 0.3
    assert data["maf"] == 0.25


def test_process_variant_snp(mock_variant_snp):
    """测试处理SNP变异位点"""
    result = va.process_variant(mock_variant_snp)

    assert result is not None
    assert result.chrom == "chr1"
    assert result.pos == 1000
    assert result.alleles == "A,T"
    assert result.missing == 0.2  # 1/5
    assert result.het == 0.5  # 2/4 (排除缺失)
    assert result.maf == 0.4


def test_process_variant_non_snp(mock_variant_non_snp):
    """测试处理非SNP变异位点"""
    result = va.process_variant(mock_variant_non_snp)
    assert result is None


def test_process_variant_multi_alt(mock_variant_multi_alt):
    """测试处理多等位基因变异位点"""
    result = va.process_variant(mock_variant_multi_alt)
    assert result is None


def test_process_variant_batch(mock_variants):
    """测试批处理变异位点"""
    results = va.process_variant_batch(mock_variants)

    # 应该只处理3个变异位点
    assert len(results) == 3

    # 检查第一个结果
    assert results[0].chrom == "chr1"
    assert results[0].pos == 1000
    assert results[0].alleles == "A,T"

    # 检查第二个结果
    assert results[1].chrom == "chr1"
    assert results[1].pos == 2000
    assert results[1].alleles == "A,TG"

    # 检查第三个结果
    assert results[2].chrom == "chr2"
    assert results[2].pos == 1500
    assert results[2].alleles == "C,G"


def test_write_results():
    """测试写入结果"""
    # 创建模拟CSV写入器和结果
    output = io.StringIO()
    writer = csv.DictWriter(
        output,
        fieldnames=["chrom", "pos", "alleles", "missing", "het", "maf"],
        delimiter="\t",
    )
    writer.writeheader()

    results = [
        va.VcfStats("chr1", 1000, "A,T", 0.2, 0.5, 0.4),
        va.VcfStats("chr2", 1500, "C,G", 0.4, 0.33, 0.25),
    ]

    # 调用函数
    va.write_results(writer, results)

    # 检查输出
    output_lines = output.getvalue().strip().split("\n")
    assert len(output_lines) == 3  # 标题行 + 2个结果
    assert "chr1\t1000\tA,T\t0.2\t0.5\t0.4" in output_lines[1]
    assert "chr2\t1500\tC,G\t0.4\t0.33\t0.25" in output_lines[2]


def test_process_completed_futures():
    """测试处理已完成的Future"""
    # 创建模拟Future和写入器
    future1 = MagicMock()
    future1.done.return_value = True
    future1.result.return_value = [va.VcfStats("chr1", 1000, "A,T", 0.2, 0.5, 0.4)]

    future2 = MagicMock()
    future2.done.return_value = False

    future3 = MagicMock()
    future3.done.return_value = True
    future3.result.return_value = [va.VcfStats("chr2", 1500, "C,G", 0.4, 0.33, 0.25)]

    futures = {future1: True, future2: True, future3: True}

    # 创建模拟写入器
    writer = MagicMock()

    # 调用函数
    processed = va.process_completed_futures(futures, writer)

    # 检查结果
    assert processed == 2  # 应该处理了2个Future
    assert len(futures) == 1  # 应该剩下1个Future
    assert future2 in futures  # 未完成的Future应该保留
    assert writer.writerow.call_count == 2  # 应该写入了2个结果


@patch("panel.simple_vcf_stats_large.open_vcf_file")
@patch("panel.simple_vcf_stats_large.process_variants_in_parallel")
def test_stream_process_vcf(mock_process, mock_open_vcf):
    """测试流式处理VCF文件"""
    # 设置模拟对象
    mock_vcf = MagicMock()
    mock_open_vcf.return_value = mock_vcf
    mock_process.return_value = 100  # 模拟处理了100个变异位点

    # 创建临时文件
    with tempfile.NamedTemporaryFile(suffix=".tsv") as tmp:
        output_file = Path(tmp.name)

        # 调用函数
        variant_count, elapsed_time = va.stream_process_vcf(
            Path("dummy.vcf"), output_file, num_processes=2, batch_size=50
        )

        # 检查结果
        assert variant_count == 100
        assert elapsed_time > 0

        # 验证调用
        mock_open_vcf.assert_called_once()
        mock_process.assert_called_once()


@patch("panel.simple_vcf_stats_large.ProcessPoolExecutor")
def test_process_variants_in_parallel(mock_executor):
    """测试并行处理变异位点"""
    # 设置模拟对象
    mock_vcf = MagicMock()
    mock_vcf.__iter__.return_value = [
        MockVariant("chr1", 1000, "A", ["T"], [0, 1, 1, 2, 3], 0.4),
        MockVariant("chr2", 1500, "C", ["G"], [0, 3, 1, 2, 3], 0.25),
    ]

    mock_writer = MagicMock()
    mock_progress = MagicMock()

    # 设置executor
    executor_instance = MagicMock()
    mock_executor.return_value.__enter__.return_value = executor_instance

    # 设置future
    future = MagicMock()
    future.done.return_value = True
    future.result.return_value = [va.VcfStats("chr1", 1000, "A,T", 0.2, 0.5, 0.4)]
    executor_instance.submit.return_value = future

    # 调用函数
    count = va.process_variants_in_parallel(mock_vcf, mock_writer, 2, 1, mock_progress)

    # 检查结果
    assert count == 2
    assert mock_progress.call_count == 2
    assert executor_instance.submit.call_count == 2
