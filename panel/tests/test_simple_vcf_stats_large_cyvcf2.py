#!/usr/bin/env python3
"""
VCF分析工具的测试模块

此模块包含对vcf_analyzer.py中函数的单元测试。
"""

import csv
import io
import os
import tempfile
from pathlib import Path
from unittest.mock import ANY, MagicMock, patch

import numpy as np
import pytest

# 导入被测试的模块
import panel.simple_vcf_stats_large_cyvcf2 as va


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
def mock_variant_data_snp():
    """创建模拟SNP变异位点数据"""
    return va.VariantData(
        chrom="chr1",
        pos=1000,
        ref="A",
        alt=["T"],
        gt_types=[0, 1, 1, 2, 3],
        aaf=0.4,
        is_snp=True,
    )


@pytest.fixture
def mock_variant_data_non_snp():
    """创建模拟非SNP变异位点数据"""
    return va.VariantData(
        chrom="chr1",
        pos=2000,
        ref="A",
        alt=["TG"],
        gt_types=[0, 1, 2, 0, 1],
        aaf=0.3,
        is_snp=False,
    )


@pytest.fixture
def mock_variant_data_multi_alt():
    """创建模拟多等位基因变异位点数据"""
    return va.VariantData(
        chrom="chr1",
        pos=3000,
        ref="A",
        alt=["T", "G"],
        gt_types=[0, 1, 2, 0, 1],
        aaf=0.3,
        is_snp=True,
    )


@pytest.fixture
def mock_variant_data_list():
    """创建模拟变异位点数据列表"""
    return [
        # 正常SNP
        va.VariantData(
            chrom="chr1",
            pos=1000,
            ref="A",
            alt=["T"],
            gt_types=[0, 1, 1, 2, 3],
            aaf=0.4,
            is_snp=True,
        ),
        # 非SNP (将被跳过)
        va.VariantData(
            chrom="chr1",
            pos=2000,
            ref="A",
            alt=["TG"],
            gt_types=[0, 1, 2, 0, 1],
            aaf=0.3,
            is_snp=False,
        ),
        # 多个替代等位基因 (将被跳过)
        va.VariantData(
            chrom="chr1",
            pos=3000,
            ref="A",
            alt=["T", "G"],
            gt_types=[0, 1, 2, 0, 1],
            aaf=0.3,
            is_snp=True,
        ),
        # 正常SNP，有缺失值
        va.VariantData(
            chrom="chr2",
            pos=1500,
            ref="C",
            alt=["G"],
            gt_types=[0, 3, 1, 2, 3],
            aaf=0.25,
            is_snp=True,
        ),
    ]


def test_variant_to_data():
    """测试变异位点转换为数据"""
    variant = MockVariant("chr1", 1000, "A", ["T"], [0, 1, 1, 2, 3], 0.4)
    data = va.variant_to_data(variant)

    assert data.chrom == "chr1"
    assert data.pos == 1000
    assert data.ref == "A"
    assert data.alt == ["T"]
    assert data.gt_types == [0, 1, 1, 2, 3]
    assert data.aaf == 0.4
    assert data.is_snp == True


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


def test_process_variant_snp(mock_variant_data_snp):
    """测试处理SNP变异位点"""
    result = va.process_variant(mock_variant_data_snp)

    assert result is not None
    assert result.chrom == "chr1"
    assert result.pos == 1000
    assert result.alleles == "A,T"
    assert result.missing == 0.2  # 1/5
    assert result.het == 0.5  # 2/4 (排除缺失)
    assert result.maf == 0.4


def test_process_variant_non_snp(mock_variant_data_non_snp):
    """测试处理非SNP变异位点"""
    result = va.process_variant(mock_variant_data_non_snp)
    assert result is None


def test_process_variant_multi_alt(mock_variant_data_multi_alt):
    """测试处理多等位基因变异位点"""
    result = va.process_variant(mock_variant_data_multi_alt)
    assert result is None


def test_process_variant_batch(mock_variant_data_list):
    """测试批处理变异位点"""
    results = va.process_variant_batch(mock_variant_data_list)

    # 应该只处理2个变异位点（第1个和第4个）
    assert len(results) == 2

    # 检查第一个结果
    assert results[0].chrom == "chr1"
    assert results[0].pos == 1000
    assert results[0].alleles == "A,T"

    # 检查第二个结果
    assert results[1].chrom == "chr2"
    assert results[1].pos == 1500
    assert results[1].alleles == "C,G"


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


@patch("panel.simple_vcf_stats_large_cyvcf2.open_vcf_file")
@patch("panel.simple_vcf_stats_large_cyvcf2.process_variants_in_parallel")
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
            Path("dummy.vcf"),
            output_file,
            num_processes=2,
            batch_size=50,
            show_progress=False,
        )

        # 检查结果
        assert variant_count == 100
        assert elapsed_time > 0

        # 验证调用
        mock_open_vcf.assert_called_once()
        mock_process.assert_called_once_with(mock_vcf, ANY, 2, 50, False)


@patch("panel.simple_vcf_stats_large_cyvcf2.ProcessPoolExecutor")
@patch("panel.simple_vcf_stats_large_cyvcf2.print_progress")
@patch("panel.simple_vcf_stats_large_cyvcf2.variant_to_data")
def test_process_variants_in_parallel(
    mock_variant_to_data, mock_print_progress, mock_executor
):
    """测试并行处理变异位点"""
    # 设置模拟对象
    mock_vcf = MagicMock()
    mock_vcf.__iter__.return_value = [
        MockVariant("chr1", 1000, "A", ["T"], [0, 1, 1, 2, 3], 0.4),
        MockVariant("chr2", 1500, "C", ["G"], [0, 3, 1, 2, 3], 0.25),
    ]

    # 设置variant_to_data的返回值
    mock_variant_to_data.side_effect = [
        va.VariantData("chr1", 1000, "A", ["T"], [0, 1, 1, 2, 3], 0.4, True),
        va.VariantData("chr2", 1500, "C", ["G"], [0, 3, 1, 2, 3], 0.25, True),
    ]

    mock_writer = MagicMock()

    # 设置executor
    executor_instance = MagicMock()
    mock_executor.return_value.__enter__.return_value = executor_instance

    # 设置future
    future = MagicMock()
    future.done.return_value = True
    future.result.return_value = [va.VcfStats("chr1", 1000, "A,T", 0.2, 0.5, 0.4)]
    executor_instance.submit.return_value = future

    # 调用函数
    count = va.process_variants_in_parallel(mock_vcf, mock_writer, 2, 1, True)

    # 检查结果
    assert count == 2
    assert mock_print_progress.call_count == 2
    assert executor_instance.submit.call_count == 2
    assert mock_variant_to_data.call_count == 2


def test_print_progress(capsys):
    """测试打印进度"""
    # 调用函数
    va.print_progress(10000)

    # 捕获输出
    captured = capsys.readouterr()
    assert "已处理 10,000 个变异位点" in captured.out

    # 测试不满足间隔的情况
    va.print_progress(9999)
    captured = capsys.readouterr()
    assert captured.out == ""  # 不应该有输出


@pytest.fixture
def sample_vcf():
    """创建一个临时的VCF文件用于测试"""
    with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as tmp:
        # 写入VCF头部信息
        tmp.write(
            b"""##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3\tSample4\tSample5
chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/1\t0/1\t1/1\t0/0
chr1\t200\t.\tG\tC\t.\t.\t.\tGT\t0/1\t0/1\t1/1\t1/1\t./.
chr2\t150\t.\tT\tA,G\t.\t.\t.\tGT\t0/0\t0/1\t0/2\t1/1\t2/2
"""
        )

    yield tmp.name
    # 测试后删除临时文件
    os.unlink(tmp.name)


def test_vcf_reading(sample_vcf):
    """测试使用cyvcf2读取VCF文件"""
    try:
        import cyvcf2
    except ImportError:
        pytest.skip("cyvcf2库未安装，跳过测试")

    # 打开VCF文件
    vcf = cyvcf2.VCF(sample_vcf, gts012=True)

    # 测试样本信息
    samples = vcf.samples
    assert len(samples) == 5
    assert samples == ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5"]

    # 读取所有变异位点
    variants = list(vcf)

    # 测试变异位点数量
    assert len(variants) == 3

    # 测试第一个变异位点
    v1 = variants[0]
    assert v1.CHROM == "chr1"
    assert v1.POS == 100
    assert v1.REF == "A"
    assert v1.ALT == ["T"]
    assert v1.is_snp

    # 测试基因型
    gt_types1 = v1.gt_types  # 0=HOM_REF, 1=HET, 2=HOM_ALT, 3=UNKNOWN
    expected_gt_types1 = np.array([0, 1, 1, 2, 0])
    assert np.array_equal(gt_types1, expected_gt_types1)

    # 测试第二个变异位点（包含缺失值）
    v2 = variants[1]
    assert v2.CHROM == "chr1"
    assert v2.POS == 200
    assert v2.REF == "G"
    assert v2.ALT == ["C"]
    assert v2.is_snp

    # 测试基因型（包含缺失值）
    gt_types2 = v2.gt_types
    expected_gt_types2 = np.array([1, 1, 2, 2, 3])  # 最后一个是缺失值
    assert np.array_equal(gt_types2, expected_gt_types2)

    # 测试第三个变异位点（多等位基因）
    v3 = variants[2]
    assert v3.CHROM == "chr2"
    assert v3.POS == 150
    assert v3.REF == "T"
    assert v3.ALT == ["A", "G"]  # 多等位基因
    assert v3.is_snp  # 多等位基因被视为非SNP

    # 测试基因型（多等位基因）
    gt_types3 = v3.gt_types
    expected_gt_types3 = np.array(
        [0, 1, 1, 2, 2]
    )  # cyvcf2将多等位基因的基因型映射为0,1,2
    assert np.array_equal(gt_types3, expected_gt_types3)

    # 测试等位基因频率
    assert 0 <= v1.aaf <= 1
    assert 0 <= v2.aaf <= 1
    assert 0 <= v3.aaf <= 1

    print("VCF文件读取测试通过！")
