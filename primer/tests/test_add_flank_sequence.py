#!/usr/bin/env python3
"""
Tests for add_flank_sequence.py
"""

import os
import sys
import tempfile
from pathlib import Path

import pandas as pd
import pytest

# Add parent directory to path to import the module
sys.path.insert(0, str(Path(__file__).parent.parent))

from add_flank_sequence import (
    add_probe_sequences,
    extract_probe_sequence,
    replace_position_with_alleles,
)


class TestExtractProbeSequence:
    """Test extract_probe_sequence function"""

    def test_extract_probe_sequence_success(self):
        """Test successful sequence extraction"""
        # Create a temporary FASTA file
        fasta_content = ">chr1\nATCGATCGATCGATCGATCG"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
            f.write(fasta_content)
            fasta_file = f.name

        try:
            # Test extraction - pysam uses 0-based coordinates
            sequence = extract_probe_sequence(fasta_file, "chr1", 5, 15)
            assert sequence == "TCGATCGATC"  # positions 5-15 (0-based)
        finally:
            os.unlink(fasta_file)

    def test_extract_probe_sequence_invalid_chrom(self):
        """Test handling of invalid chromosome"""
        # Create a temporary FASTA file
        fasta_content = ">chr1\nATCGATCGATCGATCGATCG"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
            f.write(fasta_content)
            fasta_file = f.name

        try:
            # Test with invalid chromosome - should return N's
            sequence = extract_probe_sequence(fasta_file, "chr2", 5, 15)
            assert sequence == "N" * 10  # Should return N's for invalid chrom
        finally:
            os.unlink(fasta_file)


class TestReplacePositionWithAlleles:
    """Test replace_position_with_alleles function"""

    def test_replace_position_valid_alleles(self):
        """Test replacing position with valid alleles"""
        sequence = "ATCGATCGAT"
        result = replace_position_with_alleles(sequence, 3, "A/T")
        assert result == "ATC[A/T]ATCGAT"

    def test_replace_position_invalid_alleles(self):
        """Test handling of invalid alleles format"""
        sequence = "ATCGATCGAT"
        result = replace_position_with_alleles(sequence, 3, "A")
        assert result == sequence  # Should return original sequence

    def test_replace_position_out_of_bounds(self):
        """Test handling of out-of-bounds position"""
        sequence = "ATCGATCGAT"
        result = replace_position_with_alleles(sequence, 15, "A/T")
        assert result == sequence  # Should return original sequence

    def test_replace_position_edge_cases(self):
        """Test edge cases for position replacement"""
        sequence = "ATCGATCGAT"

        # Test first position
        result = replace_position_with_alleles(sequence, 0, "A/T")
        assert result == "[A/T]TCGATCGAT"

        # Test last position (position 9 in 10-character string)
        result = replace_position_with_alleles(sequence, 9, "A/T")
        assert result == "ATCGATCGA[A/T]"  # Corrected expectation


class TestAddProbeSequences:
    """Test add_probe_sequences function"""

    def test_add_probe_sequences_success(self):
        """Test successful addition of probe sequences"""
        # Create sample data
        data = {
            "chrom": ["chr1", "chr1"],
            "pos": [1000, 2000],
            "probe_start": [995, 1995],
            "probe_end": [1005, 2005],
            "alleles": ["A/T", "G/C"],
        }
        df = pd.DataFrame(data)

        # Create temporary files
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            df.to_csv(f.name, sep="\t", index=False)
            table_file = f.name

        fasta_content = ">chr1\n" + "N" * 3000  # Create a long sequence
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
            f.write(fasta_content)
            fasta_file = f.name

        with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as f:
            output_file = f.name

        try:
            # Run the function
            add_probe_sequences(table_file, fasta_file, output_file)

            # Check output
            result_df = pd.read_csv(output_file, sep="\t")
            assert "probe_sequence" in result_df.columns
            assert len(result_df) == 2

            # Check that sequences were added
            assert not result_df["probe_sequence"].isna().any()

        finally:
            # Clean up
            os.unlink(table_file)
            os.unlink(fasta_file)
            os.unlink(output_file)

    def test_add_probe_sequences_missing_columns(self):
        """Test handling of missing required columns"""
        # Create data with missing columns
        data = {"chrom": ["chr1"], "pos": [1000], "alleles": ["A/T"]}
        df = pd.DataFrame(data)

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            df.to_csv(f.name, sep="\t", index=False)
            table_file = f.name

        fasta_content = ">chr1\nATCGATCGAT"
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
            f.write(fasta_content)
            fasta_file = f.name

        with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as f:
            output_file = f.name

        try:
            # Should raise ValueError for missing columns
            with pytest.raises(ValueError, match="Missing required columns"):
                add_probe_sequences(table_file, fasta_file, output_file)
        finally:
            # Clean up
            os.unlink(table_file)
            os.unlink(fasta_file)
            os.unlink(output_file)

    def test_add_probe_sequences_csv_format(self):
        """Test handling of CSV format input"""
        # Create sample data
        data = {
            "chrom": ["chr1"],
            "pos": [1000],
            "probe_start": [995],
            "probe_end": [1005],
            "alleles": ["A/T"],
        }
        df = pd.DataFrame(data)

        # Create temporary files
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            df.to_csv(f.name, index=False)
            table_file = f.name

        fasta_content = ">chr1\n" + "N" * 2000
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
            f.write(fasta_content)
            fasta_file = f.name

        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
            output_file = f.name

        try:
            # Run the function
            add_probe_sequences(table_file, fasta_file, output_file)

            # Check output
            result_df = pd.read_csv(output_file)
            assert "probe_sequence" in result_df.columns
            assert len(result_df) == 1

        finally:
            # Clean up
            os.unlink(table_file)
            os.unlink(fasta_file)
            os.unlink(output_file)


class TestIntegration:
    """Integration tests"""

    def test_full_workflow(self):
        """Test the complete workflow with realistic data"""
        # Create realistic test data
        data = {
            "chrom": ["chr1", "chr1", "chr2"],
            "pos": [1000, 2000, 1500],
            "probe_start": [995, 1995, 1495],
            "probe_end": [1005, 2005, 1505],
            "alleles": ["A/T", "G/C", "T/A"],
        }
        df = pd.DataFrame(data)

        # Create a realistic FASTA file
        fasta_content = ">chr1\n" + "ATCG" * 1000 + "\n>chr2\n" + "GCTA" * 1000

        # Create temporary files
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            df.to_csv(f.name, sep="\t", index=False)
            table_file = f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
            f.write(fasta_content)
            fasta_file = f.name

        with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as f:
            output_file = f.name

        try:
            # Run the function
            add_probe_sequences(table_file, fasta_file, output_file)

            # Check output
            result_df = pd.read_csv(output_file, sep="\t")

            # Verify structure
            assert "probe_sequence" in result_df.columns
            assert len(result_df) == 3

            # Verify that all sequences contain the allele format
            for seq in result_df["probe_sequence"]:
                assert "[" in seq and "]" in seq
                assert "/" in seq

            # Verify original data is preserved
            for col in ["chrom", "pos", "probe_start", "probe_end", "alleles"]:
                assert col in result_df.columns

        finally:
            # Clean up
            os.unlink(table_file)
            os.unlink(fasta_file)
            os.unlink(output_file)


if __name__ == "__main__":
    pytest.main([__file__])
