import typer
import csv
from Bio import SeqIO
import pandas as pd

from pathlib import Path

gff_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

def main(genome: Path, gff: Path, chrom_map: Path) -> None:
    genome_list = []
    chrom_map_dict = {}
    with open(chrom_map) as map_inf:
        for eachline in map_inf:
            eachline_list = eachline.strip().split()
            chrom_map_dict[eachline_list[0]] = eachline_list[1]
    for record in SeqIO.parse(genome, format='fasta'):
        record.id = chrom_map_dict[record.id]
        genome_list.append(record)

    out_genome_name = f'rename.{genome.name}'
    out_genome_path = genome.parent / out_genome_name    

    SeqIO.write(genome_list, out_genome_path, 'fasta')

    gff_df = pd.read_table(gff, header=None, names=gff_columns, comment='#')
    chrom_map_df = pd.read_table(chrom_map, header=None, names=['seqname', 'new_seqname'])
    gff_df = gff_df.merge(chrom_map_df)
    out_gff_name = f'rename.{gff.name}'
    out_gff_path = gff.parent / out_gff_name
    gff_df.to_csv(out_gff_path, sep='\t', header=False, columns=['new_seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'], index=False, quoting=csv.QUOTE_NONE)


if __name__ == "__main__":
    typer.run(main)
