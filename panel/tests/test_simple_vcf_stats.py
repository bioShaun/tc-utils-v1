# test_vcf_stats.py
import os
import tempfile
from pathlib import Path

import pandas as pd
import pysam
import pytest

from panel.simple_vcf_stats import (
    AlleleStats,
    calculate_allele_stats,
    calculate_rates,
    process_variant,
    process_vcf,
)


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


def test_calculate_allele_stats():
    """测试calculate_allele_stats函数"""

    # 测试用例1: 正常基因型
    samples1 = [
        {"GT": (0, 0)},  # 纯合参考
        {"GT": (0, 1)},  # 杂合
        {"GT": (1, 1)},  # 纯合替代
    ]
    stats1 = calculate_allele_stats(samples1)
    assert stats1.missing_count == 0
    assert stats1.het_count == 1
    assert stats1.ref_count == 3
    assert stats1.alt_count == 3

    # 测试用例2: 包含缺失值
    samples2 = [
        {"GT": (0, 0)},
        {"GT": (None, None)},  # 完全缺失
        {"GT": (0, 1)},
    ]
    stats2 = calculate_allele_stats(samples2)
    assert stats2.missing_count == 1
    assert stats2.het_count == 1
    assert stats2.ref_count == 3
    assert stats2.alt_count == 1

    # 测试用例3: 部分缺失值
    samples3 = [
        {"GT": (0, None)},  # 部分缺失
        {"GT": (1, 1)},
    ]
    stats3 = calculate_allele_stats(samples3)
    assert stats3.missing_count == 0
    assert stats3.het_count == 0
    assert stats3.ref_count == 1
    assert stats3.alt_count == 2


def test_calculate_rates():
    """测试calculate_rates函数"""
    # 测试用例1: 正常情况
    stats1 = AlleleStats(missing_count=1, het_count=2, ref_count=6, alt_count=4)
    rates1 = calculate_rates(stats1, 10)
    assert rates1.missing_rate == 0.1  # 1/10
    assert rates1.het_rate == pytest.approx(0.222, abs=0.001)  # 2/9
    assert rates1.maf == 0.4  # 4/10

    # 测试用例2: 全部缺失
    stats2 = AlleleStats(missing_count=5, het_count=0, ref_count=0, alt_count=0)
    rates2 = calculate_rates(stats2, 5)
    assert rates2.missing_rate == 1.0
    assert rates2.het_rate == 0.0
    assert rates2.maf == 0.0

    # 测试用例3: 全部杂合
    stats3 = AlleleStats(missing_count=0, het_count=5, ref_count=5, alt_count=5)
    rates3 = calculate_rates(stats3, 5)
    assert rates3.missing_rate == 0.0
    assert rates3.het_rate == 1.0
    assert rates3.maf == 0.5


def test_process_variant(sample_vcf):
    """测试process_variant函数"""
    vcf = pysam.VariantFile(sample_vcf)
    records = list(vcf)

    # 测试双等位基因位点
    result1 = process_variant(records[0])
    assert result1 is not None
    assert result1["chrom"] == "chr1"
    assert result1["pos"] == 100
    assert result1["alleles"] == "A,T"
    assert result1["missing"] == 0.0
    assert result1["het"] == 0.4
    assert result1["maf"] == 0.4

    # 测试部分缺失的位点
    result2 = process_variant(records[1])
    assert result2 is not None
    assert result2["chrom"] == "chr1"
    assert result2["pos"] == 200
    assert result2["missing"] == 0.2
    assert result2["het"] == 0.5
    assert pytest.approx(result2["maf"]) == 0.25

    # 测试多等位基因位点
    result3 = process_variant(records[2])
    assert result3 is None  # 应该返回None，因为不是双等位基因位点


def test_process_vcf(sample_vcf):
    """测试process_vcf函数"""
    # 创建临时输出文件
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp_out:
        output_file = tmp_out.name

    # 处理VCF文件
    df = process_vcf(Path(sample_vcf), Path(output_file))

    # 验证结果
    assert len(df) == 2  # 只有2个双等位基因位点
    assert list(df.columns) == [
        "chrom",
        "pos",
        "alleles",
        "missing",
        "het",
        "maf",
    ]  # 验证列名

    # 验证第一个位点的计算结果
    row1 = df.iloc[0]
    assert row1["chrom"] == "chr1"
    assert row1["pos"] == 100
    assert row1["alleles"] == "A,T"
    assert row1["missing"] == 0.0
    assert row1["het"] == 0.4
    assert row1["maf"] == 0.4

    # 验证第二个位点的计算结果
    row2 = df.iloc[1]
    assert row2["chrom"] == "chr1"
    assert row2["pos"] == 200
    assert row2["alleles"] == "G,C"
    assert row2["missing"] == 0.2
    assert row2["het"] == 0.5
    assert pytest.approx(row2["maf"]) == 0.25

    # 验证输出文件是否存在
    assert os.path.exists(output_file)

    # 验证输出文件内容
    df_from_file = pd.read_csv(output_file)
    assert len(df_from_file) == len(df)

    # 清理
    os.unlink(output_file)
