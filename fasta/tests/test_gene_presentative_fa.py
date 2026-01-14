
import pytest
import typer
from pathlib import Path
from Bio import SeqIO
from fasta.gene_presentative_fa_v2 import gene_presentative_fa

# ... fixtures ...

def test_invalid_fasta_format(tmp_path):
    """Test error when regex fails and no map provided."""
    fasta = tmp_path / "bad.fa"
    fasta.write_text(">TR1 no_gene_info\nATGC\n")
    
    with pytest.raises(typer.Exit):
        gene_presentative_fa(fasta)

def test_missing_input(tmp_path):
    """Test validation for missing files."""
    with pytest.raises(typer.Exit):
        gene_presentative_fa(tmp_path / "nonexistent.fa")
