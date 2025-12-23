import os
import tempfile
from pathlib import Path

import pytest

from draft.vcf_genotype_stats import get_genotype_class, process_vcf


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
    os.unlink(tmp.name)


class TestGetGenotypeClass:
    """测试 get_genotype_class 函数"""

    def test_hom_ref(self):
        assert get_genotype_class(0, 0) == "hom_ref"

    def test_hom_alt(self):
        assert get_genotype_class(1, 1) == "hom_alt"
        assert get_genotype_class(2, 2) == "hom_alt"

    def test_het(self):
        assert get_genotype_class(0, 1) == "het"
        assert get_genotype_class(1, 2) == "het"
        assert get_genotype_class(0, 2) == "het"

    def test_missing(self):
        assert get_genotype_class(-1, -1) == "missing"
        assert get_genotype_class(-1, 0) == "missing"
        assert get_genotype_class(0, -1) == "missing"
        assert get_genotype_class(-1, 1) == "missing"


def test_process_vcf_basic(sample_vcf):
    """测试 process_vcf 函数基本功能"""
    result = process_vcf(Path(sample_vcf))

    # 验证总体统计结果
    # chr1:100 - 2 hom_ref (Sample1, Sample5)
    # chr1:200 - 0 hom_ref
    # chr2:150 - 1 hom_ref (Sample1)
    assert result["hom_ref"] == 3
    # chr1:100 - 2 het, chr1:200 - 2 het, chr2:150 - 2 het (0/1, 0/2)
    assert result["het"] == 6
    # chr1:100 - 1 hom_alt, chr1:200 - 2 hom_alt, chr2:150 - 2 hom_alt (1/1, 2/2)
    assert result["hom_alt"] == 5
    assert result["missing"] == 1  # 1*(./.) at chr1:200


def test_process_vcf_with_output(sample_vcf):
    """测试 process_vcf 函数输出CSV文件"""
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp_out:
        output_file = tmp_out.name

    try:
        process_vcf(Path(sample_vcf), Path(output_file))

        # 验证输出文件存在
        assert os.path.exists(output_file)

        # 读取CSV文件验证内容
        import pandas as pd

        df = pd.read_csv(output_file)
        assert len(df) == 3  # 3个位点

        # 验证列名包含 het_rate 和 missing_rate
        assert "het_rate" in df.columns
        assert "missing_rate" in df.columns

        # 验证基本列存在
        assert "CHROM" in df.columns
        assert "POS" in df.columns
        assert "hom_ref" in df.columns
        assert "het" in df.columns
        assert "hom_alt" in df.columns
        assert "missing" in df.columns

    finally:
        if os.path.exists(output_file):
            os.unlink(output_file)


def test_het_rate_calculation(sample_vcf):
    """测试杂合率计算"""
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp_out:
        output_file = tmp_out.name

    try:
        process_vcf(Path(sample_vcf), Path(output_file))

        import pandas as pd

        df = pd.read_csv(output_file)

        # chr1:100 - 5 samples: 2 het, 2 hom_ref, 1 hom_alt, 0 missing
        # het_rate = 2/5 = 40.0, missing_rate = 0%
        row1 = df[df["POS"] == 100].iloc[0]
        assert row1["het"] == 2
        assert row1["het_rate"] == 40.0
        assert row1["missing_rate"] == 0.0

        # chr1:200 - 5 samples: 2 het, 0 hom_ref, 2 hom_alt, 1 missing
        # het_rate = 2/4 = 50.0, missing_rate = 1/5 = 20.0
        row2 = df[df["POS"] == 200].iloc[0]
        assert row2["het"] == 2
        assert row2["het_rate"] == 50.0
        assert row2["missing_rate"] == 20.0

    finally:
        if os.path.exists(output_file):
            os.unlink(output_file)


def test_missing_rate_calculation(sample_vcf):
    """测试缺失率计算"""
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp_out:
        output_file = tmp_out.name

    try:
        process_vcf(Path(sample_vcf), Path(output_file))

        import pandas as pd

        df = pd.read_csv(output_file)

        # chr1:200 - 5 samples: 1 missing
        # missing_rate = 1/5 = 20.0
        row = df[df["POS"] == 200].iloc[0]
        assert row["missing"] == 1
        assert row["missing_rate"] == 20.0

    finally:
        if os.path.exists(output_file):
            os.unlink(output_file)


def test_all_missing_genotypes():
    """测试全部为缺失基因型的情况"""
    with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as tmp:
        tmp.write(
            b"""##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2
chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t./.\t./.
"""
        )
        vcf_path = tmp.name

    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp_out:
        output_file = tmp_out.name

    try:
        result = process_vcf(Path(vcf_path), Path(output_file))

        assert result["missing"] == 2
        assert result["hom_ref"] == 0
        assert result["het"] == 0
        assert result["hom_alt"] == 0

        import pandas as pd

        df = pd.read_csv(output_file)
        assert df.iloc[0]["missing_rate"] == 100.0
        assert df.iloc[0]["het_rate"] == 0.0

    finally:
        os.unlink(vcf_path)
        if os.path.exists(output_file):
            os.unlink(output_file)


def test_all_het_genotypes():
    """测试全部为杂合基因型的情况"""
    with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as tmp:
        tmp.write(
            b"""##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3
chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\t0/1
"""
        )
        vcf_path = tmp.name

    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp_out:
        output_file = tmp_out.name

    try:
        result = process_vcf(Path(vcf_path), Path(output_file))

        assert result["het"] == 3
        assert result["missing"] == 0

        import pandas as pd

        df = pd.read_csv(output_file)
        assert df.iloc[0]["het_rate"] == 100.0
        assert df.iloc[0]["missing_rate"] == 0.0

    finally:
        os.unlink(vcf_path)
        if os.path.exists(output_file):
            os.unlink(output_file)
