"""
Pytest configuration for primer tests
"""

import sys
from pathlib import Path

import pytest

# Add the primer directory to the path for imports
primer_dir = Path(__file__).parent.parent
sys.path.insert(0, str(primer_dir))


@pytest.fixture
def sample_fasta_file(tmp_path):
    """Create a sample FASTA file for testing"""
    fasta_file = tmp_path / "test.fa"
    fasta_content = ">chr1\nATCGATCGATCGATCGATCG\n>chr2\nGCTAGCTAGCTAGCTAGCTA"
    fasta_file.write_text(fasta_content)
    return str(fasta_file)


@pytest.fixture
def sample_table_data():
    """Sample table data for testing"""
    return {
        "chrom": ["chr1", "chr1", "chr2"],
        "pos": [1000, 2000, 1500],
        "probe_start": [995, 1995, 1495],
        "probe_end": [1005, 2005, 1505],
        "alleles": ["A/T", "G/C", "T/A"],
    }


@pytest.fixture
def sample_table_file(tmp_path, sample_table_data):
    """Create a sample table file for testing"""
    import pandas as pd

    table_file = tmp_path / "test.tsv"
    df = pd.DataFrame(sample_table_data)
    df.to_csv(table_file, sep="\t", index=False)
    return str(table_file)
