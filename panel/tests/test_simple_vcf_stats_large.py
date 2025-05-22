#!/usr/bin/env python3
"""
VCF分析工具的测试模块

此模块包含对vcf_analyzer.py中函数的单元测试，适用于pysam版本。
"""

import csv
import io
import os
import tempfile
from pathlib import Path
from unittest.mock import ANY, MagicMock, mock_open, patch

import numpy as np
import pytest

# 导入被测试的模块
import panel.simple_vcf_stats_large as va


class MockVariantRecord:
    """模拟pysam.VariantRecord对象"""

    def __init__(self, chrom, pos, ref, alts, samples_data):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alts = alts
        self.samples = {}

        for sample_name, gt in samples_data.items():
            self.samples[sample_name] = {"GT": gt}


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
    # 创建模拟pysam变异位点
    samples_data = {
        "Sample1": (0, 0),
        "Sample2": (0, 1),
        "Sample3": (0, 1),
        "Sample4": (1, 1),
        "Sample5": (None, None),
    }
    variant = MockVariantRecord("chr1", 1000, "A", ("T",), samples_data)
    samples = ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5"]

    # 调用函数
    data = va.variant_to_data(variant, samples)

    # 验证结果
    assert data.chrom == "chr1"
    assert data.pos == 1000
    assert data.ref == "A"
    assert data.alt == ["T"]
    assert data.gt_types == [0, 1, 1, 2, 3]  # 0=HOM_REF, 1=HET, 2=HOM_ALT, 3=UNKNOWN
    assert 0.3 <= data.aaf <= 0.5  # 大约应该是0.4，但允许有一定误差
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
        mock_process.assert_called_once_with(mock_vcf, ANY, ANY, 2, 50, False)


@patch("panel.simple_vcf_stats_large.ProcessPoolExecutor")
@patch("panel.simple_vcf_stats_large.print_progress")
@patch("panel.simple_vcf_stats_large.variant_to_data")
def test_process_variants_in_parallel(
    mock_variant_to_data, mock_print_progress, mock_executor
):
    """测试并行处理变异位点"""
    # 设置模拟对象
    mock_vcf = MagicMock()
    mock_vcf.header.samples = ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5"]

    # 创建模拟变异位点
    samples_data1 = {
        "Sample1": (0, 0),
        "Sample2": (0, 1),
        "Sample3": (0, 1),
        "Sample4": (1, 1),
        "Sample5": (0, 0),
    }
    variant1 = MockVariantRecord("chr1", 1000, "A", ("T",), samples_data1)

    samples_data2 = {
        "Sample1": (0, 1),
        "Sample2": (0, 1),
        "Sample3": (1, 1),
        "Sample4": (1, 1),
        "Sample5": (None, None),
    }
    variant2 = MockVariantRecord("chr1", 2000, "G", ("C",), samples_data2)

    mock_vcf.__iter__.return_value = [variant1, variant2]

    # 设置variant_to_data的返回值
    mock_variant_to_data.side_effect = [
        va.VariantData("chr1", 1000, "A", ["T"], [0, 1, 1, 2, 0], 0.4, True),
        va.VariantData("chr1", 2000, "G", ["C"], [1, 1, 2, 2, 3], 0.5, True),
    ]

    mock_writer = MagicMock()
    samples = ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5"]

    # 设置executor
    executor_instance = MagicMock()
    mock_executor.return_value.__enter__.return_value = executor_instance

    # 设置future
    future = MagicMock()
    future.done.return_value = True
    future.result.return_value = [va.VcfStats("chr1", 1000, "A,T", 0.2, 0.5, 0.4)]
    executor_instance.submit.return_value = future

    # 调用函数
    count = va.process_variants_in_parallel(mock_vcf, samples, mock_writer, 2, 1, True)

    # 检查结果
    assert count == 2
    assert mock_print_progress.call_count == 2
    assert executor_instance.submit.call_count == 2
    assert mock_variant_to_data.call_count == 2

    # 验证variant_to_data调用参数
    mock_variant_to_data.assert_any_call(variant1, samples)
    mock_variant_to_data.assert_any_call(variant2, samples)


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


def test_vcf_reading_pysam(sample_vcf):
    """测试使用pysam读取VCF文件"""
    try:
        import pysam
    except ImportError:
        pytest.skip("pysam库未安装，跳过测试")

    # 打开VCF文件
    vcf = pysam.VariantFile(sample_vcf)

    # 测试样本信息
    samples = list(vcf.header.samples)
    assert len(samples) == 5
    assert samples == ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5"]

    # 读取所有变异位点
    variants = list(vcf)

    # 测试变异位点数量
    assert len(variants) == 3

    # 测试第一个变异位点
    v1 = variants[0]
    assert v1.chrom == "chr1"
    assert v1.pos == 100
    assert v1.ref == "A"
    assert v1.alts == ("T",)

    # 测试基因型
    gt_types1 = []
    for sample in samples:
        gt = v1.samples[sample]["GT"]
        if gt == (0, 0):
            gt_types1.append(0)  # HOM_REF
        elif gt == (0, 1) or gt == (1, 0):
            gt_types1.append(1)  # HET
        elif gt == (1, 1):
            gt_types1.append(2)  # HOM_ALT
        else:
            gt_types1.append(3)  # UNKNOWN

    expected_gt_types1 = [0, 1, 1, 2, 0]
    assert gt_types1 == expected_gt_types1

    # 测试第二个变异位点（包含缺失值）
    v2 = variants[1]
    assert v2.chrom == "chr1"
    assert v2.pos == 200
    assert v2.ref == "G"
    assert v2.alts == ("C",)

    # 测试基因型（包含缺失值）
    gt_types2 = []
    for sample in samples:
        if "GT" not in v2.samples[sample] or None in v2.samples[sample]["GT"]:
            gt_types2.append(3)  # UNKNOWN
        else:
            gt = v2.samples[sample]["GT"]
            if gt == (0, 0):
                gt_types2.append(0)  # HOM_REF
            elif gt == (0, 1) or gt == (1, 0):
                gt_types2.append(1)  # HET
            elif gt == (1, 1):
                gt_types2.append(2)  # HOM_ALT
            else:
                gt_types2.append(3)  # UNKNOWN

    expected_gt_types2 = [1, 1, 2, 2, 3]  # 最后一个是缺失值
    assert gt_types2 == expected_gt_types2

    # 测试第三个变异位点（多等位基因）
    v3 = variants[2]
    assert v3.chrom == "chr2"
    assert v3.pos == 150
    assert v3.ref == "T"
    assert v3.alts == ("A", "G")  # 多等位基因

    # 测试基因型（多等位基因）
    gt_types3 = []
    for sample in samples:
        gt = v3.samples[sample]["GT"]
        if gt == (0, 0):
            gt_types3.append(0)  # HOM_REF
        elif gt[0] != gt[1]:
            gt_types3.append(1)  # HET
        elif gt[0] > 0 and gt[1] > 0:
            gt_types3.append(2)  # HOM_ALT
        else:
            gt_types3.append(3)  # UNKNOWN

    expected_gt_types3 = [0, 1, 1, 2, 2]
    assert gt_types3 == expected_gt_types3

    print("VCF文件读取测试通过！")
