#!/usr/bin/env python3
"""Tests for emsdb/add_db_id.py"""

import tempfile
from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from emsdb.add_db_id import (
    app,
    merge_tables,
    read_table_a_lazy,
    read_table_b_lazy,
)

runner = CliRunner()


# Sample data for Table A (with header)
TABLE_A_CONTENT = """id\tvariant\tchrom\tpos\trefer\talt\ttype\timpact\tgene
125084465\tChr1A_590095386\t1A\t590095386\tT\tG\t['missense_variant']\tMODERATE\tChr1A.g10998
124548800\tChr1A_25356524\t1A\t25356524\tC\tT\t['missense_variant']\tMODERATE\tChr1A.g01171
124554275\tChr1A_27585075\t1A\t27585075\tC\tT\t['synonymous_variant']\tLOW\tChr1A.g01262
"""

# Sample data for Table B (without header)
TABLE_B_CONTENT = """1A\t590095386\tT\tG\tmissense_variant\tMODERATE\tChr1A.g10998\tChr1A.g10998.m1\t2/2\tc.1332T>G\tp.Asp444Glu
1A\t25356524\tC\tT\tmissense_variant\tMODERATE\tChr1A.g01171\tChr1A.g01171.m1\t14/25\tc.2648G>A\tp.Gly883Asp
1A\t27585075\tC\tT\tsynonymous_variant\tLOW\tChr1A.g01262\tChr1A.g01262.m1\t3/5\tc.387C>T\tp.Phe129Phe
1A\t999999\tA\tG\tmissense_variant\tMODERATE\tChr1A.g99999\tChr1A.g99999.m1\t1/1\tc.100A>G\tp.Met1Val
"""

# Sample data for Table B (with header)
TABLE_B_WITH_HEADER = """chrom\tpos\trefer\talt\ttype\timpact\tgene\ttranscript\texon_rank\tcds_pos\tprotein_pos
1A\t590095386\tT\tG\tmissense_variant\tMODERATE\tChr1A.g10998\tChr1A.g10998.m1\t2/2\tc.1332T>G\tp.Asp444Glu
1A\t25356524\tC\tT\tmissense_variant\tMODERATE\tChr1A.g01171\tChr1A.g01171.m1\t14/25\tc.2648G>A\tp.Gly883Asp
"""


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def table_a_file(temp_dir):
    """Create a temporary Table A file."""
    path = temp_dir / "table_a.tsv"
    path.write_text(TABLE_A_CONTENT)
    return path


@pytest.fixture
def table_b_file(temp_dir):
    """Create a temporary Table B file (no header)."""
    path = temp_dir / "table_b.tsv"
    path.write_text(TABLE_B_CONTENT)
    return path


@pytest.fixture
def table_b_with_header_file(temp_dir):
    """Create a temporary Table B file (with header)."""
    path = temp_dir / "table_b_header.tsv"
    path.write_text(TABLE_B_WITH_HEADER)
    return path


class TestReadTableALazy:
    """Tests for read_table_a_lazy function."""

    def test_read_table_a_basic(self, table_a_file):
        """Test basic reading of Table A."""
        lf = read_table_a_lazy(str(table_a_file), "\t")
        df = lf.collect()

        assert df.shape[0] == 3
        assert set(df.columns) == {"id", "chrom", "pos", "refer", "alt"}

    def test_read_table_a_columns(self, table_a_file):
        """Test that only required columns are selected."""
        lf = read_table_a_lazy(str(table_a_file), "\t")
        df = lf.collect()

        # Should not contain other columns from original file
        assert "variant" not in df.columns
        assert "type" not in df.columns
        assert "impact" not in df.columns

    def test_read_table_a_pos_type(self, table_a_file):
        """Test that pos column is cast to Int64."""
        lf = read_table_a_lazy(str(table_a_file), "\t")
        df = lf.collect()

        assert df["pos"].dtype == pl.Int64


class TestReadTableBLazy:
    """Tests for read_table_b_lazy function."""

    def test_read_table_b_no_header(self, table_b_file):
        """Test reading Table B without header."""
        lf = read_table_b_lazy(str(table_b_file), "\t", has_header=False)
        df = lf.collect()

        assert df.shape[0] == 4
        expected_cols = {
            "chrom",
            "pos",
            "refer",
            "alt",
            "type",
            "impact",
            "gene",
            "transcript",
            "exon_rank",
            "cds_pos",
            "protein_pos",
        }
        assert set(df.columns) == expected_cols

    def test_read_table_b_with_header(self, table_b_with_header_file):
        """Test reading Table B with header."""
        lf = read_table_b_lazy(str(table_b_with_header_file), "\t", has_header=True)
        df = lf.collect()

        assert df.shape[0] == 2
        assert "chrom" in df.columns
        assert "transcript" in df.columns

    def test_read_table_b_pos_type(self, table_b_file):
        """Test that pos column is cast to Int64."""
        lf = read_table_b_lazy(str(table_b_file), "\t", has_header=False)
        df = lf.collect()

        assert df["pos"].dtype == pl.Int64


class TestMergeTables:
    """Tests for merge_tables function."""

    def test_merge_basic(self, table_a_file, table_b_file):
        """Test basic merge of two tables."""
        lf_a = read_table_a_lazy(str(table_a_file), "\t")
        lf_b = read_table_b_lazy(str(table_b_file), "\t", has_header=False)

        result = merge_tables(lf_a, lf_b).collect()

        # Should match 3 records (one in table_b has no match in table_a)
        assert result.shape[0] == 3

    def test_merge_output_columns(self, table_a_file, table_b_file):
        """Test that merge produces correct output columns."""
        lf_a = read_table_a_lazy(str(table_a_file), "\t")
        lf_b = read_table_b_lazy(str(table_b_file), "\t", has_header=False)

        result = merge_tables(lf_a, lf_b).collect()

        expected_cols = [
            "id",
            "type",
            "impact",
            "gene",
            "transcript",
            "exon_rank",
            "cds_pos",
            "protein_pos",
        ]
        assert result.columns == expected_cols

    def test_merge_id_from_table_a(self, table_a_file, table_b_file):
        """Test that id comes from Table A."""
        lf_a = read_table_a_lazy(str(table_a_file), "\t")
        lf_b = read_table_b_lazy(str(table_b_file), "\t", has_header=False)

        result = merge_tables(lf_a, lf_b).collect()

        # Check that IDs from table A are present
        ids = result["id"].to_list()
        assert 125084465 in ids
        assert 124548800 in ids
        assert 124554275 in ids

    def test_merge_annotations_from_table_b(self, table_a_file, table_b_file):
        """Test that annotations come from Table B."""
        lf_a = read_table_a_lazy(str(table_a_file), "\t")
        lf_b = read_table_b_lazy(str(table_b_file), "\t", has_header=False)

        result = merge_tables(lf_a, lf_b).collect()

        # Filter by known id and check annotation
        row = result.filter(pl.col("id") == 125084465)
        assert row["transcript"][0] == "Chr1A.g10998.m1"
        assert row["cds_pos"][0] == "c.1332T>G"
        assert row["protein_pos"][0] == "p.Asp444Glu"


class TestCLI:
    """Tests for CLI interface."""

    def test_cli_basic(self, table_a_file, table_b_file, temp_dir):
        """Test basic CLI usage."""
        output_file = temp_dir / "output.tsv"

        result = runner.invoke(
            app,
            [
                "-a",
                str(table_a_file),
                "-b",
                str(table_b_file),
                "-o",
                str(output_file),
            ],
        )

        assert result.exit_code == 0
        assert output_file.exists()

    def test_cli_with_header(self, table_a_file, table_b_with_header_file, temp_dir):
        """Test CLI with --table-b-header flag."""
        output_file = temp_dir / "output.tsv"

        result = runner.invoke(
            app,
            [
                "-a",
                str(table_a_file),
                "-b",
                str(table_b_with_header_file),
                "-o",
                str(output_file),
                "--table-b-header",
            ],
        )

        assert result.exit_code == 0
        assert output_file.exists()

    def test_cli_output_content(self, table_a_file, table_b_file, temp_dir):
        """Test that CLI output file has correct content."""
        output_file = temp_dir / "output.tsv"

        runner.invoke(
            app,
            [
                "-a",
                str(table_a_file),
                "-b",
                str(table_b_file),
                "-o",
                str(output_file),
            ],
        )

        # Read output and verify
        df = pl.read_csv(output_file, separator="\t")
        assert df.shape[0] == 3
        assert "id" in df.columns
        assert "transcript" in df.columns

    def test_cli_missing_table_a(self, table_b_file, temp_dir):
        """Test CLI error when Table A is missing."""
        output_file = temp_dir / "output.tsv"

        result = runner.invoke(
            app,
            [
                "-a",
                "/nonexistent/file.tsv",
                "-b",
                str(table_b_file),
                "-o",
                str(output_file),
            ],
        )

        assert result.exit_code == 1

    def test_cli_missing_table_b(self, table_a_file, temp_dir):
        """Test CLI error when Table B is missing."""
        output_file = temp_dir / "output.tsv"

        result = runner.invoke(
            app,
            [
                "-a",
                str(table_a_file),
                "-b",
                "/nonexistent/file.tsv",
                "-o",
                str(output_file),
            ],
        )

        assert result.exit_code == 1

    def test_cli_help(self):
        """Test CLI help message."""
        result = runner.invoke(app, ["--help"])

        assert result.exit_code == 0
        assert "table-a" in result.output
        assert "table-b" in result.output
        assert "output" in result.output
