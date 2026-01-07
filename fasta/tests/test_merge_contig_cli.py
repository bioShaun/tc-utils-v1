from pathlib import Path

import pytest
from typer.testing import CliRunner

from fasta.merge_contig_v2 import app

runner = CliRunner()


@pytest.fixture
def sample_files(tmp_path):
    """Create sample input files for CLI tests."""
    # Create genome FASTA
    genome_fa = tmp_path / "genome.fa"
    genome_fa.write_text(">ctg1\nAAAA\n>ctg2\nTTTT\n>chr1\nGGGG\n")

    # Create contig list
    contig_list = tmp_path / "contig_list.txt"
    contig_list.write_text("ctg1\nctg2\n")

    # Create GTF
    gtf_file = tmp_path / "genes.gtf"
    gtf_file.write_text('ctg1\ttest\tgene\t1\t2\t.\t+\t.\tgene_id "G1";\n')

    # Create offset file (for gtf command)
    offset_file = tmp_path / "offset.txt"
    offset_file.write_text("contig_id\toffset\nctg1\t0\n")

    return {
        "genome_fa": genome_fa,
        "contig_list": contig_list,
        "gtf_file": gtf_file,
        "offset_file": offset_file,
    }


def test_cli_merge_help():
    """Test that --help displays help message."""
    result = runner.invoke(app, ["merge", "--help"])
    assert result.exit_code == 0
    assert "Merge genome contigs" in result.stdout


def test_cli_merge_basic(sample_files):
    """Test basic merge command execution."""
    genome = str(sample_files["genome_fa"])
    contig_list = str(sample_files["contig_list"])

    result = runner.invoke(app, ["merge", genome, contig_list])

    assert result.exit_code == 0
    assert "Creating merged contig" in result.stdout
    assert Path(genome).with_suffix(".merge_ctg.fa").exists()


def test_cli_merge_with_gtf(sample_files):
    """Test merge command with GTF update."""
    genome = str(sample_files["genome_fa"])
    contig_list = str(sample_files["contig_list"])
    gtf = str(sample_files["gtf_file"])

    result = runner.invoke(app, ["merge", genome, contig_list, "--gtf", gtf])

    assert result.exit_code == 0
    assert "GTF file provided" in result.stdout
    assert Path(gtf).with_suffix(".merge_ctg.gtf").exists()


def test_cli_merge_file_not_found(tmp_path):
    """Test proper error exit when file is missing."""
    result = runner.invoke(app, ["merge", "nonexistent.fa", "list.txt"])

    assert result.exit_code == 1
    assert "file not found" in result.stdout


def test_cli_gtf_basic(sample_files):
    """Test gtf command execution."""
    gtf = str(sample_files["gtf_file"])
    offset = str(sample_files["offset_file"])

    result = runner.invoke(app, ["gtf", gtf, offset])

    assert result.exit_code == 0
    assert "Written merged GTF" in result.stdout
    assert Path(gtf).with_suffix(".merge_ctg.gtf").exists()
