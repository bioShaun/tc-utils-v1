import os
import tempfile
from pathlib import Path

import pandas as pd
import pytest
from rich import print

from chip.genetic_inheritance import VCFGeneticAnalyzer

# 构造一个最小的VCF内容（2个位点，3个样本）
VCF_CONTENT = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tF1\tM1\tC1
1\t100\trs1\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t0/1
1\t200\trs2\tC\tT\t60\tPASS\t.\tGT\t0/1\t1/1\t0/1
"""

# 构造一个最小的家系组合表
FAMILY_CONTENT = "father,mother,child\nF1,M1,C1\n"


@pytest.fixture
def temp_vcf_file():
    with tempfile.NamedTemporaryFile(delete=False, suffix=".vcf") as f:
        f.write(VCF_CONTENT.encode())
        temp_path = Path(f.name)
    yield temp_path
    os.remove(temp_path)


@pytest.fixture
def temp_family_file():
    with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as f:
        f.write(FAMILY_CONTENT.encode())
        temp_path = Path(f.name)
    yield temp_path
    os.remove(temp_path)


def test_vcfgeneticanalyzer_basic(temp_vcf_file, temp_family_file):
    analyzer = VCFGeneticAnalyzer()
    vcf_df = analyzer.load_vcf_data(temp_vcf_file)
    assert not vcf_df.empty
    fam_df = analyzer.load_family_combinations(temp_family_file)
    assert not fam_df.empty
    results = analyzer.analyze_family_consistency()
    assert not results.empty
    assert "consistency_rate" in results.columns
    # 检查一致性比例在合理范围
    assert (
        results["consistency_rate"].iloc[0] >= 0.0
        and results["consistency_rate"].iloc[0] <= 1.0
    )
    # 检查详细报告和保存结果不会报错
    analyzer.generate_detailed_report(output_dir=".")
    with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as out:
        analyzer.save_results(Path(out.name))
    # 清理
    os.remove(out.name)


def make_vcf_content(father_gt, mother_gt, child_gt):
    return f"""##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tF1\tM1\tC1
1\t100\trs1\tA\tG\t50\tPASS\t.\tGT\t{father_gt}\t{mother_gt}\t{child_gt}
"""


def run_family_case(father_gt, mother_gt, child_gt, expected_consistent):
    # 生成临时VCF和家系文件
    vcf_content = make_vcf_content(father_gt, mother_gt, child_gt)
    with tempfile.NamedTemporaryFile(delete=False, suffix=".vcf") as f:
        f.write(vcf_content.encode())
        vcf_path = Path(f.name)
    with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as f:
        f.write(FAMILY_CONTENT.encode())
        fam_path = Path(f.name)
    analyzer = VCFGeneticAnalyzer()
    vcf_df = analyzer.load_vcf_data(vcf_path)
    fam_df = analyzer.load_family_combinations(fam_path)
    results = analyzer.analyze_family_consistency()
    assert not results.empty
    rate = results["consistency_rate"].iloc[0]
    if expected_consistent:
        assert rate == 1.0, f"{father_gt},{mother_gt},{child_gt} 应该一致"
    else:
        assert rate == 0.0, f"{father_gt},{mother_gt},{child_gt} 应该不一致"
    os.remove(vcf_path)
    os.remove(fam_path)


def test_all_family_scenarios():
    # 所有二倍体基因型组合
    gts = ["0/0", "0/1", "1/1"]
    # 合法子代组合
    valid = {
        ("0/0", "0/0"): ["0/0"],
        ("0/0", "0/1"): ["0/0", "0/1"],
        ("0/0", "1/1"): ["0/1"],
        ("0/1", "0/0"): ["0/0", "0/1"],
        ("0/1", "0/1"): ["0/0", "0/1", "1/1"],
        ("0/1", "1/1"): ["0/1", "1/1"],
        ("1/1", "0/0"): ["0/1"],
        ("1/1", "0/1"): ["0/1", "1/1"],
        ("1/1", "1/1"): ["1/1"],
    }
    for fg in gts:
        for mg in gts:
            for cg in gts:
                expected = cg in valid[(fg, mg)]
                run_family_case(fg, mg, cg, expected)
