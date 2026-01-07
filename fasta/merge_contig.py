from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import fire
from pandas import DataFrame
from pathlib import PurePath, Path
import pandas as pd


def merge_contig_gtf(gtf_file, ctg_offset, new_name='chrUn'):
    gtf_file = Path(gtf_file)
    file_sfx = gtf_file.suffix
    new_gtf_file = gtf_file.with_suffix(f'.merge_ctg{file_sfx}')
    ctg_offset_df = pd.read_table(ctg_offset, index_col=0)
    ctg_offset_df.index = [str(each) for each in ctg_offset_df.index]
    new_gtf_inf = open(new_gtf_file, 'w')
    with open(gtf_file) as gtf_inf:
        for eachline in gtf_inf:
            if eachline.startswith('#') or eachline.strip() == '':
                new_gtf_inf.write(eachline)
                continue
            eachline_inf = eachline.strip().split('\t')
            chrom = eachline_inf[0]
            try:
                start = int(eachline_inf[3])
            except IndexError:
                print(eachline)
            end = int(eachline_inf[4])
            if chrom in ctg_offset_df.index:
                start = start + ctg_offset_df.loc[chrom, 'offset']
                end = end + ctg_offset_df.loc[chrom, 'offset']
                chrom = new_name
            eachline_inf[0] = chrom
            eachline_inf[3] = str(start)
            eachline_inf[4] = str(end)
            output_line_str = '\t'.join(eachline_inf)
            new_gtf_inf.write(f'{output_line_str}\n')
    new_gtf_inf.close()


def merge_contig_fa(genome_fa, congtig_list, n_sep=100,
                    merge_name='chrUn'):
    '''
    merge contigs to a super contig in genome file.
    each contig is seperated with a number of N in super contig.
    output a new genome file with the super contig and
    a table with each contig's offset in super contig
    '''
    contig_df = pd.read_table(congtig_list, index_col=0, header=None)
    contig_df.index = [str(each) for each in contig_df.index]
    n_sep_str = 'N' * n_sep
    contig_seq_list = []
    congtig_offset_dict = {}
    offset = 0
    all_seq_list = []
    for seq_record in SeqIO.parse(genome_fa, "fasta"):
        if seq_record.id in contig_df.index:
            contig_seq_list.append(str(seq_record.seq))
            congtig_offset_dict.setdefault(
                'contig_id', []).append(seq_record.id)
            congtig_offset_dict.setdefault(
                'offset', []).append(offset)
            offset += len(seq_record.seq) + n_sep
        else:
            all_seq_list.append(seq_record)
    merged_contig_seq = Seq(n_sep_str.join(contig_seq_list))
    merged_contig_seq_rd = SeqRecord(id=merge_name,
                                     seq=merged_contig_seq,
                                     description='')
    all_seq_list.append(merged_contig_seq_rd)
    genome_fa = PurePath(genome_fa)
    genome_merge_ctg_fa = genome_fa.with_suffix('.merge_ctg.fa')
    SeqIO.write(all_seq_list, genome_merge_ctg_fa, "fasta")
    ctg_offset_file = genome_fa.with_suffix('.ctg.offset.txt')
    ctg_offset_df = DataFrame(congtig_offset_dict)
    ctg_offset_df.to_csv(ctg_offset_file, sep='\t', index=False,
                         columns=['contig_id', 'offset'])
    return genome_merge_ctg_fa, ctg_offset_file


def merge_contig(genome_fa, contig_list, n_sep=100,
                 merge_name='chrUn', gtf_file=None):
    # merge contig fastas
    _, ctg_offset_file = merge_contig_fa(genome_fa, contig_list, n_sep=100,
                                         merge_name=merge_name)
    # change contig names in gtf
    if gtf_file is not None:
        merge_contig_gtf(gtf_file, ctg_offset_file, merge_name)


if __name__ == '__main__':
    fire.Fire({'merge_contig': merge_contig,
               'merge_contig_gtf': merge_contig_gtf})
