import re
import sys
from pathlib import PurePath

import fire
import gtfparse
import pandas as pd
from Bio import Seq, SeqIO


def gene_presentative_fa(gffread_fa, gtf=None, gene_tr_map=None):
    if gtf is not None:
        gtf_df = gtfparse.read_gtf(gtf)
        tr_gtf_df = gtf_df[gtf_df.feature == "exon"].loc[
            :, ["transcript_id", "gene_id"]
        ]
        tr_gtf_df = tr_gtf_df.drop_duplicates().set_index("transcript_id")
    elif gene_tr_map is not None:
        tr_gtf_df = pd.read_csv(
            gene_tr_map, header=None, index_col=1, names=["gene_id"], sep="\t"
        )
    else:
        tr_gtf_df = None
    gene_fa_dict = dict()
    for seq_record in SeqIO.parse(gffread_fa, "fasta"):
        # replace dot in seq, diamond mkindex fail
        seq_record.seq = Seq.Seq(str(seq_record.seq).replace(".", "*"))
        if tr_gtf_df is None:
            if re.search(r"gene=(\S+)", seq_record.description):
                gene_id = re.search(r"gene=(\S+)", seq_record.description).groups()[0]
            else:
                sys.exit(
                    "Wrong fasta format. [>transcript_id gene=gene_id]: {}".format(
                        seq_record.description
                    )
                )
        else:
            if seq_record.id in tr_gtf_df.index:
                gene_id = str(tr_gtf_df.loc[seq_record.id].gene_id)
            else:
                continue
        seq_len = len(seq_record.seq)
        seq_record.id = gene_id
        seq_record.description = ""
        if gene_id in gene_fa_dict:
            if gene_fa_dict[gene_id][0] >= seq_len:
                continue
        gene_fa_dict[gene_id] = [seq_len, seq_record]

    seq_record_list = []
    for gene_i in gene_fa_dict:
        seq_record_list.append(gene_fa_dict[gene_i][1])
    gffread_presentative_fa = PurePath(gffread_fa).with_suffix(".gene.fa")
    SeqIO.write(seq_record_list, gffread_presentative_fa, "fasta")


if __name__ == "__main__":
    fire.Fire(gene_presentative_fa)
