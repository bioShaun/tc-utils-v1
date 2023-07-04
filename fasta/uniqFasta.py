import typer

from Bio import SeqIO
from pathlib import Path


def uniqFasta(fa_file: Path, out_file: Path) -> None:
    seqIdSet = set()
    seqList = []
    for seq_record in SeqIO.parse(fa_file, "fasta"):
        if seq_record.id in seqIdSet:
            continue
        seqIdSet.add(seq_record.id)
        seqList.append(seq_record)
    SeqIO.write(seqList, out_file, "fasta")


if __name__ == "__main__":
    typer.run(uniqFasta)
