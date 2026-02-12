#!/usr/bin/env python3
"""Tests for primer/vcf2primer.py."""

from pathlib import Path

import pandas as pd
import pytest
import typer
from vcf2primer import main


def write_fasta(tmp_path: Path) -> Path:
    """Create a test FASTA file."""
    fasta_path = tmp_path / "ref.fa"
    fasta_path.write_text(">chr1\nACGTACGTAC\n", encoding="utf-8")
    return fasta_path


def write_vcf(tmp_path: Path, records: list[tuple[str, int, str, str, str]]) -> Path:
    """Create a test VCF file with provided records."""
    vcf_path = tmp_path / "input.vcf"
    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    lines.extend(
        f"{chrom}\t{pos}\t{site_id}\t{ref}\t{alt}\t.\t.\t."
        for chrom, pos, site_id, ref, alt in records
    )
    vcf_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return vcf_path


def write_id_file(tmp_path: Path, ids: list[str], file_name: str = "ids.txt") -> Path:
    """Create an ID list file."""
    id_path = tmp_path / file_name
    id_path.write_text("\n".join(ids) + "\n", encoding="utf-8")
    return id_path


def test_main_without_id_uses_chrom_pos_when_all_vcf_ids_are_dot(tmp_path: Path) -> None:
    """When all VCF IDs are '.', output names should use chrom_pos."""
    ref = write_fasta(tmp_path)
    vcf = write_vcf(
        tmp_path,
        [
            ("chr1", 3, ".", "G", "A"),
            ("chr1", 6, ".", "C", "T"),
        ],
    )
    out = tmp_path / "out.tsv"

    main(vcf, ref, out, flank_size=2)

    result = pd.read_table(out)
    assert list(result["name"]) == ["chr1_3", "chr1_6"]
    assert list(result["sequence"]) == ["AC[G/A]TA", "TA[C/T]GT"]


def test_main_filters_by_id_file_when_all_vcf_ids_are_dot(tmp_path: Path) -> None:
    """ID filtering should work with chrom_pos IDs when all VCF IDs are '.'."""
    ref = write_fasta(tmp_path)
    vcf = write_vcf(
        tmp_path,
        [
            ("chr1", 3, ".", "G", "A"),
            ("chr1", 6, ".", "C", "T"),
        ],
    )
    id_file = write_id_file(tmp_path, ["chr1_6"])
    out = tmp_path / "out.tsv"

    main(vcf, ref, out, flank_size=2, id_file=id_file)

    result = pd.read_table(out)
    assert list(result["name"]) == ["chr1_6"]
    assert len(result) == 1


def test_main_filters_by_original_id_when_vcf_ids_are_not_all_dot(tmp_path: Path) -> None:
    """ID filtering should use original VCF IDs when they are not all '.'."""
    ref = write_fasta(tmp_path)
    vcf = write_vcf(
        tmp_path,
        [
            ("chr1", 3, "rs1", "G", "A"),
            ("chr1", 6, ".", "C", "T"),
        ],
    )
    id_file = write_id_file(tmp_path, ["rs1"])
    out = tmp_path / "out.tsv"

    main(vcf, ref, out, flank_size=2, id_file=id_file)

    result = pd.read_table(out)
    assert list(result["name"]) == ["rs1"]
    assert len(result) == 1


def test_main_does_not_fallback_to_chrom_pos_for_mixed_ids(tmp_path: Path) -> None:
    """For mixed IDs, chrom_pos should not be used as filter key."""
    ref = write_fasta(tmp_path)
    vcf = write_vcf(
        tmp_path,
        [
            ("chr1", 3, "rs1", "G", "A"),
            ("chr1", 6, ".", "C", "T"),
        ],
    )
    id_file = write_id_file(tmp_path, ["chr1_3"])
    out = tmp_path / "out.tsv"

    with pytest.raises(typer.Exit) as exc_info:
        main(vcf, ref, out, flank_size=2, id_file=id_file)

    assert exc_info.value.exit_code == 1


def test_main_ignores_blank_lines_in_id_file(tmp_path: Path) -> None:
    """Blank lines in the ID file should be ignored."""
    ref = write_fasta(tmp_path)
    vcf = write_vcf(
        tmp_path,
        [
            ("chr1", 3, ".", "G", "A"),
            ("chr1", 6, ".", "C", "T"),
        ],
    )
    id_file = write_id_file(tmp_path, ["", "chr1_3", ""])
    out = tmp_path / "out.tsv"

    main(vcf, ref, out, flank_size=2, id_file=id_file)

    result = pd.read_table(out)
    assert list(result["name"]) == ["chr1_3"]


def test_main_raises_when_id_file_not_found(tmp_path: Path) -> None:
    """Missing ID file should fail with exit code 1."""
    ref = write_fasta(tmp_path)
    vcf = write_vcf(tmp_path, [("chr1", 3, ".", "G", "A")])
    out = tmp_path / "out.tsv"
    missing_id_file = tmp_path / "missing_ids.txt"

    with pytest.raises(typer.Exit) as exc_info:
        main(vcf, ref, out, flank_size=2, id_file=missing_id_file)

    assert exc_info.value.exit_code == 1


def test_main_extracts_insertion_and_deletion_sequences_correctly(tmp_path: Path) -> None:
    """Insertion/deletion sites should produce correct flanking sequences."""
    ref = write_fasta(tmp_path)
    vcf = write_vcf(
        tmp_path,
        [
            ("chr1", 4, ".", "T", "TA"),   # insertion
            ("chr1", 3, ".", "GTA", "G"),  # deletion
        ],
    )
    out = tmp_path / "out.tsv"

    main(vcf, ref, out, flank_size=2)

    result = pd.read_table(out)
    sequence_by_name = dict(zip(result["name"], result["sequence"], strict=True))
    assert sequence_by_name["chr1_4"] == "CG[T/TA]AC"
    assert sequence_by_name["chr1_3"] == "AC[GTA/G]CG"


def test_main_filters_indel_by_id_file(tmp_path: Path) -> None:
    """ID filtering should work for indel records."""
    ref = write_fasta(tmp_path)
    vcf = write_vcf(
        tmp_path,
        [
            ("chr1", 4, ".", "T", "TA"),
            ("chr1", 3, ".", "GTA", "G"),
        ],
    )
    id_file = write_id_file(tmp_path, ["chr1_4"])
    out = tmp_path / "out.tsv"

    main(vcf, ref, out, flank_size=2, id_file=id_file)

    result = pd.read_table(out)
    assert list(result["name"]) == ["chr1_4"]
    assert list(result["sequence"]) == ["CG[T/TA]AC"]
