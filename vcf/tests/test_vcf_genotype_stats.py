import csv
import pytest
from pathlib import Path
from typer.testing import CliRunner
from vcf.vcf_genotype_stats import app, get_genotype_class, process_vcf, print_summary

@pytest.mark.parametrize("a1, a2, expected", [
    (0, 0, "hom_ref"),
    (1, 1, "hom_alt"),
    (2, 2, "hom_alt"),
    (0, 1, "het"),
    (1, 0, "het"),
    (1, 2, "het"),
    (-1, -1, "missing"),
    (-1, 0, "missing"),
    (0, -1, "missing"),
    (1, -1, "missing"),
    (-1, 1, "missing"),
])
def test_get_genotype_class(a1, a2, expected):
    """Test the genotype classification logic."""
    assert get_genotype_class(a1, a2) == expected

@pytest.fixture
def sample_vcf(tmp_path):
    """Create a small sample VCF file."""
    vcf_path = tmp_path / "sample.vcf"
    content = (
        "##fileformat=VCFv4.2\n"
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
        "##contig=<ID=chr1,length=1000>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3\n"
        "chr1\t100\t.\tA\tG\t100\tPASS\t.\tGT\t0/0\t0/1\t1/1\n"
        "chr1\t200\t.\tT\tC\t100\tPASS\t.\tGT\t0/0\t./.\t1/1\n"
        "chr1\t300\t.\tC\tG\t100\tPASS\t.\tGT\t0/1\t0/1\t0/1\n"
    )
    vcf_path.write_text(content)
    return vcf_path

def test_process_vcf(sample_vcf):
    """Test the process_vcf function."""
    stats = process_vcf(sample_vcf)
    
    assert stats["hom_ref"] == 2
    assert stats["het"] == 4
    assert stats["hom_alt"] == 2
    assert stats["missing"] == 1

def test_process_vcf_with_csv(sample_vcf, tmp_path):
    """Test the process_vcf function with CSV output."""
    output_csv = tmp_path / "stats.csv"
    process_vcf(sample_vcf, output_path=output_csv)
    
    assert output_csv.exists()
    
    with open(output_csv, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        
    assert len(rows) == 3
    assert rows[0]["CHROM"] == "chr1"
    assert rows[0]["POS"] == "100"
    assert int(rows[0]["hom_ref"]) == 1
    assert int(rows[0]["het"]) == 1
    assert int(rows[0]["hom_alt"]) == 1
    assert int(rows[0]["missing"]) == 0
    
    assert rows[1]["POS"] == "200"
    assert int(rows[1]["missing"]) == 1

def test_print_summary(capsys):
    """Test the print_summary function."""
    stats = {"hom_ref": 10, "het": 5, "hom_alt": 5, "missing": 2}
    print_summary(stats)
    captured = capsys.readouterr()
    assert "Genotype Statistics Summary" in captured.out
    assert "Homozygous REF (0/0):" in captured.out
    assert "10" in captured.out

def test_analyze_cli(sample_vcf):
    """Test the analyze command via CliRunner."""
    runner = CliRunner()
    result = runner.invoke(app, [str(sample_vcf)])
    assert result.exit_code == 0
    assert "Genotype Statistics Summary" in result.stdout

def test_analyze_cli_error():
    """Test CLI error handling with non-existent file."""
    runner = CliRunner()
    result = runner.invoke(app, ["non_existent.vcf"])
    assert result.exit_code != 0

def test_import():
    """Verify that the module can be imported."""
    try:
        from vcf import vcf_genotype_stats
    except ImportError:
        assert False, "Failed to import vcf.vcf_genotype_stats"