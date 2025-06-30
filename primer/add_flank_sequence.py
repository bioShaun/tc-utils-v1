from pathlib import Path
from typing import Optional

import pandas as pd
import pysam
import typer


def extract_probe_sequence(
    fasta_file: str, chrom: str, probe_start: int, probe_end: int
) -> str:
    """
    Extract sequence from probe_start to probe_end from FASTA file.

    Args:
        fasta_file: Path to FASTA file
        chrom: Chromosome name
        probe_start: Probe start position (0-based)
        probe_end: Probe end position (0-based)

    Returns:
        Probe sequence string
    """
    with pysam.FastaFile(fasta_file) as fa:
        try:
            sequence = fa.fetch(chrom, probe_start, probe_end)
            return sequence
        except (ValueError, KeyError) as e:
            print(
                f"Warning: Could not fetch sequence for {chrom}:{probe_start}-{probe_end}: {e}"
            )
            return "N" * (probe_end - probe_start)


def replace_position_with_alleles(sequence: str, pos_in_seq: int, alleles: str) -> str:
    """
    Replace the position in sequence with alleles in [A/T] format.

    Args:
        sequence: Original sequence
        pos_in_seq: Position within the sequence (0-based)
        alleles: Alleles string like "A/T", "G/C"

    Returns:
        Sequence with position replaced by [alleles]
    """
    if pos_in_seq >= len(sequence):
        return sequence

    # Split alleles and format
    allele_list = alleles.split("/")
    if len(allele_list) != 2:
        print(f"Warning: Invalid alleles format '{alleles}', expected 'A/T' format")
        return sequence

    # Replace the position with [alleles]
    new_sequence = sequence[:pos_in_seq] + f"[{alleles}]" + sequence[pos_in_seq + 1 :]

    return new_sequence


def add_probe_sequences(table_file: str, genome_fasta: str, output_file: str) -> None:
    """
    Add probe sequences to table with genomic positions and alleles.

    Args:
        table_file: Input table file (CSV/TSV)
        genome_fasta: Genome FASTA file
        output_file: Output table file
    """
    # Read table
    if table_file.endswith(".csv"):
        df = pd.read_csv(table_file)
    else:
        df = pd.read_csv(table_file, sep="\t")

    # Check required columns
    required_cols = ["chrom", "pos", "probe_start", "probe_end", "alleles"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Add probe sequence column
    probe_sequences = []

    for idx, row in df.iterrows():
        chrom = str(row["chrom"])
        pos = int(row["pos"])
        probe_start = int(row["probe_start"])
        probe_end = int(row["probe_end"])
        alleles = str(row["alleles"])

        # Extract probe sequence
        probe_seq = extract_probe_sequence(genome_fasta, chrom, probe_start, probe_end)

        # Calculate position within probe sequence (0-based)
        # pos is 1-based, probe_start is 0-based
        pos_in_probe = pos - 1 - probe_start

        # Replace position with alleles
        modified_seq = replace_position_with_alleles(probe_seq, pos_in_probe, alleles)

        probe_sequences.append(modified_seq)

        if (idx + 1) % 100 == 0:
            print(f"Processed {idx + 1} rows...")

    # Add probe column to dataframe
    df["probe_sequence"] = probe_sequences

    # Save output
    if output_file.endswith(".csv"):
        df.to_csv(output_file, index=False)
    else:
        df.to_csv(output_file, sep="\t", index=False)

    print(f"Added probe sequences to {len(df)} rows. Output saved to {output_file}")


def main(
    table_file: str = typer.Argument(..., help="Input table file (CSV/TSV)"),
    genome_fasta: str = typer.Argument(..., help="Genome FASTA file"),
    output_file: str = typer.Argument(..., help="Output table file"),
):
    """
    Add probe sequences to table with genomic positions and alleles.

    Expected table columns:
    - chrom: Chromosome name
    - pos: 1-based position
    - probe_start: Probe start position (0-based)
    - probe_end: Probe end position (0-based)
    - alleles: Alleles in format like "A/T", "G/C"

    The script will add a 'probe_sequence' column with sequences from probe_start
    to probe_end where the position is replaced by [alleles] format.
    """
    try:
        add_probe_sequences(table_file, genome_fasta, output_file)
    except Exception as e:
        print(f"Error: {e}")
        raise typer.Exit(1)


if __name__ == "__main__":
    typer.run(main)
