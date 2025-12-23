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
    """Create a small sample VCF file with 3 samples."""
    vcf_path = tmp_path / "sample.vcf"
    # Row 1: 0/0, 0/0, 0/0 -> PASS
    # Row 2: 0/0, 0/1, 0/0 -> HET>... (1/3 = 0.33)
    # Row 3: 0/0, ./., 0/0 -> MISS>... (1/3 = 0.33)
    # Row 4: 0/1, 1/1, 1/1 -> HET>... (1/3=0.33) OR NO_REF_OR_ALT if HET threshold is high enough
    content = (
        "##fileformat=VCFv4.2\n"
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
        "##contig=<ID=chr1,length=1000>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3\n"
        "chr1\t100\t.\tA\tG\t100\tPASS\t.\tGT\t0/0\t0/0\t0/0\n"
        "chr1\t200\t.\tT\tC\t100\tPASS\t.\tGT\t0/0\t0/1\t0/0\n"
        "chr1\t300\t.\tC\tG\t100\tPASS\t.\tGT\t0/0\t./.\t0/0\n"
        "chr1\t400\t.\tG\tA\t100\tPASS\t.\tGT\t0/1\t1/1\t1/1\n"
    )
    vcf_path.write_text(content)
    return vcf_path

def test_process_vcf(sample_vcf):
    """Test the process_vcf function."""
    stats = process_vcf(sample_vcf)
    
    # 4 sites * 3 samples = 12 genotypes
    # 100: 3 Ref
    # 200: 2 Ref, 1 Het
    # 300: 2 Ref, 1 Miss
    # 400: 1 Het, 2 Alt
    # Totals: Ref=7, Het=2, Alt=2, Miss=1
    assert stats["hom_ref"] == 7
    assert stats["het"] == 2
    assert stats["hom_alt"] == 2
    assert stats["missing"] == 1

def test_process_vcf_with_csv_defaults(sample_vcf, tmp_path):
    """Test process_vcf with CSV output using default thresholds (0.1)."""
    output_csv = tmp_path / "stats_default.csv"
    process_vcf(sample_vcf, output_path=output_csv, max_missing=0.1, max_het=0.1)
    
    assert output_csv.exists()
    
    with open(output_csv, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        
    assert len(rows) == 4
    # Row 1 (100): PASS
    assert rows[0]["FILTER"] == "PASS"
    
    # Row 2 (200): Het 0.33 > 0.1 -> HET>0.1
    assert rows[1]["FILTER"] == "HET>0.1"
    
    # Row 3 (300): Miss 0.33 > 0.1 -> MISS>0.1
    assert rows[2]["FILTER"] == "MISS>0.1"
    
    # Row 4 (400): Het 0.33 > 0.1 -> HET>0.1
    assert rows[3]["FILTER"] == "HET>0.1"

def test_process_vcf_custom_thresholds(sample_vcf, tmp_path):
    """Test process_vcf with relaxed thresholds to check precedence and NO_REF_OR_ALT."""
    output_csv = tmp_path / "stats_relaxed.csv"
    # Relax thresholds to 0.4 so 0.33 doesn't trigger HET/MISS filters
    process_vcf(sample_vcf, output_path=output_csv, max_missing=0.4, max_het=0.4)
    
    with open(output_csv, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    # Row 1 (100): PASS
    assert rows[0]["FILTER"] == "PASS"
    
    # Row 2 (200): Het 0.33 <= 0.4.
    # Ref > 0 (2), Het > 0 (1), Alt = 0.
    # has_het=T, has_ref=T, has_alt=F -> NO_REF_OR_ALT
    assert rows[1]["FILTER"] == "NO_REF_OR_ALT"
    
    # Row 3 (300): Miss 0.33 <= 0.4.
    # Ref > 0 (2), Het = 0, Alt = 0. -> PASS
    assert rows[2]["FILTER"] == "PASS"
    
    # Row 4 (400): Het 0.33 <= 0.4.
    # Ref = 0, Het > 0 (1), Alt > 0 (2).
    # has_het=T, has_ref=F, has_alt=T -> NO_REF_OR_ALT
    assert rows[3]["FILTER"] == "NO_REF_OR_ALT"

def test_print_summary(capsys):
    """Test the print_summary function."""
    stats = {"hom_ref": 10, "het": 5, "hom_alt": 5, "missing": 2}
    print_summary(stats)
    captured = capsys.readouterr()
    assert "Genotype Statistics Summary" in captured.out
    assert "Homozygous REF (0/0):" in captured.out
    assert "10" in captured.out

def test_analyze_cli(sample_vcf, tmp_path):
    """Test the analyze command with new arguments."""
    runner = CliRunner()
    output_csv = tmp_path / "cli_out.csv"
    
    result = runner.invoke(app, [
        str(sample_vcf),
        "--output", str(output_csv),
        "--max-missing", "0.4",
        "--max-het", "0.4"
    ])
    
    assert result.exit_code == 0
    assert "Genotype Statistics Summary" in result.stdout
    assert output_csv.exists()
    
    # Verify content lightly to ensure args were passed
    with open(output_csv, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    # Based on test_process_vcf_custom_thresholds, row 2 (pos 200) should be NO_REF_OR_ALT
    assert rows[1]["FILTER"] == "NO_REF_OR_ALT"

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