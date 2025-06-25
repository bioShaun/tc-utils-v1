from pathlib import Path

import pandas as pd
import pytest
from typer.testing import CliRunner

from chip.genetic_inheritance import VCFGeneticAnalyzer, main

# ---- 1. 单元测试：基因型字符串解析 ----


def test_parse_vcf_genotype_cases():
    analyzer = VCFGeneticAnalyzer()
    assert analyzer.parse_vcf_genotype("0/1") == {"0", "1"}
    assert analyzer.parse_vcf_genotype("1/1") == {"1"}
    assert analyzer.parse_vcf_genotype("0|1:35:99") == {"0", "1"}
    assert analyzer.parse_vcf_genotype("1") == {"1"}
    assert analyzer.parse_vcf_genotype("./.") == set()
    assert analyzer.parse_vcf_genotype(".") == set()
    assert analyzer.parse_vcf_genotype("") == set()
    assert analyzer.parse_vcf_genotype(None) == set()


def test_resolve_alleles(tmp_path):
    analyzer = VCFGeneticAnalyzer()
    # 假设有三种等位基因
    analyzer.variant_info = pd.DataFrame(
        [
            {
                "CHROM": "1",
                "POS": 100,
                "ID": "rs1",
                "REF": "A",
                "ALT": "C,G",
                "QUAL": 50,
                "FILTER": ".",
                "INFO": ".",
                "FORMAT": "GT",
            }
        ],
        index=["1:100:A:C,G"],
    )
    assert analyzer.resolve_alleles("1:100:A:C,G", {"0", "1", "2"}) == {"A", "C", "G"}
    assert analyzer.resolve_alleles("1:100:A:C,G", {"1"}) == {"C"}
    assert analyzer.resolve_alleles("1:100:A:C,G", set()) == set()
    assert analyzer.resolve_alleles("not_exist", {"0"}) == {"0"}  # 不存在返回索引


# ---- 2. 集成测试：小VCF和家系表的全流程 ----

VCF_CONTENT = """\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tF1\tM1\tC1
1\t100\trs1\tA\tG\t100\tPASS\t.\tGT\t0/0\t1/1\t0/1
1\t200\trs2\tT\tC\t99\tPASS\t.\tGT\t0/1\t0/0\t0/1
"""


@pytest.fixture
def tiny_vcf(tmp_path):
    vcf = tmp_path / "tiny.vcf"
    with open(vcf, "w") as f:
        f.write(VCF_CONTENT)
    return vcf


@pytest.fixture
def tiny_family(tmp_path):
    fam = tmp_path / "fam.csv"
    with open(fam, "w") as f:
        f.write("father,mother,child\nF1,M1,C1\n")
    return fam


def test_load_and_analyze_full(tmp_path, tiny_vcf, tiny_family):
    analyzer = VCFGeneticAnalyzer()
    vcf_df = analyzer.load_vcf_data(tiny_vcf)
    fam_df = analyzer.load_family_combinations(tiny_family)
    assert not vcf_df.empty
    assert not fam_df.empty
    results = analyzer.analyze_family_consistency(chunk_size=1)
    assert not results.empty
    row = results.iloc[0]
    assert row["consistent_variants"] > 0
    assert row["missing_data_variants"] == 0
    # 生成报告文件
    out = tmp_path / "out.csv"
    analyzer.save_results(out)
    assert out.exists()
    df_out = pd.read_csv(out)
    assert not df_out.empty


# ---- 3. CLI测试 ----


@pytest.fixture
def cli_runner():
    return CliRunner()


def test_cli_full(tmp_path, tiny_vcf, tiny_family, cli_runner):
    out = tmp_path / "cli_out.csv"
    result = cli_runner.invoke(main, [str(tiny_vcf), str(tiny_family), str(out)])
    assert result.exit_code == 0
    assert out.exists()
    df = pd.read_csv(out)
    assert not df.empty
    # 检查命令行输出
    assert "分析完成" in result.output or "分析失败" not in result.output
