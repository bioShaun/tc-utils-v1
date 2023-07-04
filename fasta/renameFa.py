import typer

from Bio import SeqIO

from pathlib import Path


def main(fasta: Path, name_map: Path, out_fasta: Path) -> None:
    name_map_dict = dict([each.strip().split() for each in name_map.open()])
    renamed_fasta_list = []
    for record in SeqIO.parse(fasta, "fasta"):
        record.id = name_map_dict.get(record.id, record.id)
        renamed_fasta_list.append(record)
    SeqIO.write(renamed_fasta_list, out_fasta, "fasta")


if __name__ == "__main__":
    typer.run(main)
