import pytest
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from fasta.merge_contig_v2 import merge_contig_fa, merge_contig_gtf, merge_contig


@pytest.fixture
def sample_fasta(tmp_path):
    """Create a sample FASTA file with 3 contigs."""
    fasta_path = tmp_path / "genome.fa"
    records = [
        SeqRecord(Seq("AAAA"), id="ctg1", description=""),
        SeqRecord(Seq("TTTT"), id="ctg2", description=""),
        SeqRecord(Seq("CCCC"), id="ctg3", description=""),
        SeqRecord(Seq("GGGG"), id="chr1", description=""),
    ]
    SeqIO.write(records, fasta_path, "fasta")
    return fasta_path


@pytest.fixture
def sample_contig_list(tmp_path):
    """Create a contig list file for ctg1 and ctg2."""
    list_path = tmp_path / "contig_list.txt"
    list_path.write_text("ctg1\nctg2\n")
    return list_path


@pytest.fixture
def sample_gtf(tmp_path):
    """Create a sample GTF file."""
    gtf_path = tmp_path / "test.gtf"
    gtf_content = [
        "# This is a comment",
        'ctg1\ttest\tgene\t1\t2\t.\t+\t.\tgene_id "G1";',
        'ctg2\ttest\tgene\t3\t4\t.\t+\t.\tgene_id "G2";',
        'chr1\ttest\tgene\t10\t20\t.\t+\t.\tgene_id "G3";',
        "",
        'ctg1\ttest\tCDS\t3\t4\t.\t+\t.\tgene_id "G1";',
    ]
    gtf_path.write_text("\n".join(gtf_content) + "\n")
    return gtf_path


def test_merge_contig_fa(sample_fasta, sample_contig_list, tmp_path):
    """Test merging FASTA contigs."""
    genome_fa = str(sample_fasta)
    contig_list = str(sample_contig_list)

    # Run merge_contig_fa
    # Note: n_sep=100 is default, but we'll use a smaller one for easier verification
    out_fa, offset_file = merge_contig_fa(genome_fa, contig_list, n_sep=10, merge_name="chrUn")

    # Verify outputs exist
    assert Path(out_fa).exists()
    assert Path(offset_file).exists()

    # Verify FASTA content
    records = list(SeqIO.parse(out_fa, "fasta"))
    assert len(records) == 3  # chr1, ctg3, and merged chrUn (ctg1 and ctg2 were merged)

    # Check sequences
    record_dict = {r.id: str(r.seq) for r in records}
    assert "chr1" in record_dict
    assert record_dict["chr1"] == "GGGG"
    assert "ctg3" in record_dict
    assert record_dict["ctg3"] == "CCCC"
    assert "chrUn" in record_dict
    # ctg1 (AAAA) + 10 Ns + ctg2 (TTTT)
    assert record_dict["chrUn"] == "AAAA" + "N" * 10 + "TTTT"

    # Verify offset table
    offset_df = pd.read_table(offset_file)
    assert len(offset_df) == 2
    assert offset_df.iloc[0]["contig_id"] == "ctg1"
    assert offset_df.iloc[0]["offset"] == 0
    assert offset_df.iloc[1]["contig_id"] == "ctg2"
    assert offset_df.iloc[1]["offset"] == 14  # 4 (length of ctg1) + 10 (n_sep)


def test_merge_contig_gtf(sample_gtf, tmp_path):
    """Test updating GTF coordinates."""
    # Create a mock offset file
    offset_file = tmp_path / "offset.txt"
    offset_file.write_text("contig_id\toffset\nctg1\t0\nctg2\t100\n")

    # Run merge_contig_gtf
    merge_contig_gtf(str(sample_gtf), str(offset_file), new_name="chrUn")

    # Verify output file
    # suffix .gtf becomes .merge_ctg.gtf
    expected_out = sample_gtf.with_suffix(".merge_ctg.gtf")
    assert expected_out.exists()

    # Verify content
    lines = expected_out.read_text().splitlines()
    assert lines[0] == "# This is a comment"

    # ctg1 -> chrUn, start=1+0, end=2+0
    assert 'chrUn\ttest\tgene\t1\t2\t.\t+\t.\tgene_id "G1";' in lines
    # ctg2 -> chrUn, start=3+100, end=4+100
    assert 'chrUn\ttest\tgene\t103\t104\t.\t+\t.\tgene_id "G2";' in lines
    # chr1 remains unchanged
    assert 'chr1\ttest\tgene\t10\t20\t.\t+\t.\tgene_id "G3";' in lines
    # Empty line preserved
    assert "" in lines
    # Another ctg1 entry
    assert 'chrUn\ttest\tCDS\t3\t4\t.\t+\t.\tgene_id "G1";' in lines


def test_merge_contig_integration(sample_fasta, sample_contig_list, sample_gtf):
    """Test full merge_contig workflow."""
    merge_contig(
        genome_fa=str(sample_fasta),
        contig_list=str(sample_contig_list),
        n_sep=100,
        merge_name="chrUn",
        gtf_file=str(sample_gtf),
    )

    # Verify all outputs
    assert Path(sample_fasta).with_suffix(".merge_ctg.fa").exists()
    assert Path(sample_fasta).with_suffix(".ctg.offset.txt").exists()
    assert Path(sample_gtf).with_suffix(".merge_ctg.gtf").exists()


def test_merge_contig_no_gtf(sample_fasta, sample_contig_list):
    """Test merge_contig without GTF."""
    merge_contig(genome_fa=str(sample_fasta), contig_list=str(sample_contig_list), gtf_file=None)
    assert Path(sample_fasta).with_suffix(".merge_ctg.fa").exists()
    # Check that GTF output was NOT created
    assert not Path(sample_fasta).with_suffix(".merge_ctg.gtf").exists()


def test_merge_contig_fa_custom_sep(sample_fasta, sample_contig_list):
    """Test merge_contig_fa with custom separator."""
    out_fa, _ = merge_contig_fa(str(sample_fasta), str(sample_contig_list), n_sep=5)
    records = list(SeqIO.parse(out_fa, "fasta"))
    chr_un = next(r for r in records if r.id == "chrUn")
    assert str(chr_un.seq) == "AAAA" + "N" * 5 + "TTTT"


def test_merge_contig_gtf_missing_contig(sample_gtf, tmp_path):
    """Test GTF update when contig is not in offset table."""
    offset_file = tmp_path / "offset_missing.txt"
    offset_file.write_text("contig_id\toffset\nctg1\t0\n")  # ctg2 is missing

    merge_contig_gtf(str(sample_gtf), str(offset_file))
    expected_out = sample_gtf.with_suffix(".merge_ctg.gtf")
    lines = expected_out.read_text().splitlines()

    # ctg1 should be updated
    assert 'chrUn\ttest\tgene\t1\t2\t.\t+\t.\tgene_id "G1";' in lines
    # ctg2 should remain as ctg2 (since it's not in the offset table)
    assert 'ctg2\ttest\tgene\t3\t4\t.\t+\t.\tgene_id "G2";' in lines


def test_merge_contig_gtf_short_line(tmp_path):
    """Test GTF update with a short line raises IndexError with context."""
    gtf_path = tmp_path / "short.gtf"
    gtf_path.write_text("short\tline\n")

    offset_file = tmp_path / "offset.txt"
    offset_file.write_text("contig_id\toffset\nctg1\t0\n")

    with pytest.raises(IndexError, match="GTF 行列数不足"):
        merge_contig_gtf(str(gtf_path), str(offset_file))
