from pathlib import Path

import pandas as pd
from pyfaidx import Fasta

from exome import vcf2probe


def write_reference(tmp_path: Path) -> Path:
    fasta_path = tmp_path / "ref.fa"
    fasta_path.write_text(">chr1\nAACCGGTTAACC\n", encoding="utf-8")
    return fasta_path


def write_vcf(tmp_path: Path) -> Path:
    vcf_path = tmp_path / "variants.vcf"
    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "chr1\t5\t.\tG\tT\t.\tPASS\t.",
        "chr1\t7\t.\tT\tA\t.\tPASS\t.",
    ]
    vcf_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return vcf_path


def test_build_sequence_strings_returns_expected_sequences(tmp_path):
    fasta_path = write_reference(tmp_path)
    reference = Fasta(str(fasta_path))
    row = pd.Series({"chrom": "chr1", "pos": 5, "ref": "G", "alt": "T"})

    variant_seq, reference_seq = vcf2probe.build_sequence_strings(
        row, reference, half_length=2
    )

    assert variant_seq == "CC[G/T]GT"
    assert reference_seq == "CCGGT"
    reference.close()


def test_run_generates_variant_and_reference_outputs(tmp_path):
    fasta_path = write_reference(tmp_path)
    vcf_path = write_vcf(tmp_path)
    out_dir = tmp_path / "out"

    vcf2probe.run(
        vcf_path,
        fasta_path,
        out_dir,
        half_length=2,
        is_vcf=True,
        sequence_mode="both",
        show_progress=False,
    )

    table_path = out_dir / vcf2probe.TABLE_NAME
    df = pd.read_csv(table_path)
    assert set(["sequence", "reference_sequence"]).issubset(df.columns)
    first_row = df[df["pos"] == 5].iloc[0]
    assert first_row["sequence"] == "CC[G/T]GT"
    assert first_row["reference_sequence"] == "CCGGT"

    variant_fasta = (
        (out_dir / vcf2probe.VARIANT_FASTA_NAME)
        .read_text(encoding="utf-8")
        .strip()
        .splitlines()
    )
    reference_fasta = (
        (out_dir / vcf2probe.REFERENCE_FASTA_NAME)
        .read_text(encoding="utf-8")
        .strip()
        .splitlines()
    )

    assert variant_fasta[1] == "CC[G/T]GT"
    assert reference_fasta[1] == "CCGGT"


def test_reference_mode_skips_variant_fasta(tmp_path):
    fasta_path = write_reference(tmp_path)
    vcf_path = write_vcf(tmp_path)
    out_dir = tmp_path / "out_ref"

    vcf2probe.run(
        vcf_path,
        fasta_path,
        out_dir,
        half_length=2,
        is_vcf=True,
        sequence_mode="reference",
        show_progress=False,
    )

    variant_path = out_dir / vcf2probe.VARIANT_FASTA_NAME
    reference_path = out_dir / vcf2probe.REFERENCE_FASTA_NAME

    assert not variant_path.exists()
    reference_lines = reference_path.read_text(encoding="utf-8").strip().splitlines()
    assert reference_lines[1] == "CCGGT"


def test_variant_mode_skips_reference_fasta(tmp_path):
    fasta_path = write_reference(tmp_path)
    vcf_path = write_vcf(tmp_path)
    out_dir = tmp_path / "out_variant"

    vcf2probe.run(
        vcf_path,
        fasta_path,
        out_dir,
        half_length=2,
        is_vcf=True,
        sequence_mode="variant",
        show_progress=False,
    )

    reference_path = out_dir / vcf2probe.REFERENCE_FASTA_NAME
    variant_path = out_dir / vcf2probe.VARIANT_FASTA_NAME

    assert not reference_path.exists()
    variant_lines = variant_path.read_text(encoding="utf-8").strip().splitlines()
    assert variant_lines[1] == "CC[G/T]GT"
