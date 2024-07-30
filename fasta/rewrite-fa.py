from pathlib import Path

import typer
from Bio import SeqIO


def rewrite_fasta(input_file: Path, output_file: Path) -> None:
    """
    Rewrite fasta file

    Args:
        input_file (Path): Path to input fasta file.
        output_file (Path): Path to output fasta file.
    """
    complete_records = []
    for record in SeqIO.parse(input_file, "fasta"):
        if not record.seq.missing_data:
            complete_records.append(record)

    SeqIO.write(complete_records, output_file, "fasta")


if __name__ == "__main__":
    typer.run(rewrite_fasta)
