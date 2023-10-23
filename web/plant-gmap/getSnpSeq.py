from pathlib import Path
import typer
import pandas as pd
from Bio import SeqIO


COLUMN_IDX_MAP = {"varfilter": [1, 2, 3, 10]}


def get_snp_sequence(
    chrom: str,
    start: int,
    end: int,
    flank: int,
    score_type: str,
    score_file: Path,
    out_dir: Path,
    genome: Path,
    out_seq_num: int = 10,
):
    if not out_dir.is_dir():
        out_dir.mkdir()
    use_cols = COLUMN_IDX_MAP[score_type]
    df = pd.read_csv(score_file, usecols=use_cols)
    df = df[
        (df["POS"] >= start) & (df["POS"] <= end) & (df["Type"] == "SNP")
    ].drop_duplicates()
    snp_seq_list = []
    snp_seq_df_list = []
    snp_seq_out_file = (
        out_dir / f"{chrom}-{start}-{end}-flank_{flank}.snp.sequence.txt.gz"
    )
    for record in SeqIO.parse(genome, "fasta"):
        if record.id == chrom:
            for row in df.itertuples():
                ref_pos = row.POS - 1
                left_start = ref_pos - flank if ref_pos >= flank else 0
                left_seq = record.seq[left_start:ref_pos]
                right_end = ref_pos + flank + 1
                right_seq = record.seq[ref_pos + 1 : right_end]
                snp_seq = f"{left_seq}[{row.REF}/{row.ALT}]{right_seq}"
                snp_name = f"{chrom}_{row.POS}"
                snp_seq_str = f"{snp_name}\t{snp_seq}"
                if len(snp_seq_list) < out_seq_num:
                    snp_seq_list.append(snp_seq_str)
                snp_seq_df_list.append({"id": snp_name, "seq": snp_seq})
        break
    snp_seq_df = pd.DataFrame(snp_seq_df_list)
    if len(snp_seq_df) > out_seq_num:
        snp_seq_list.append("...\n请下载完整结果")
    snp_seq_df.to_csv(snp_seq_out_file, sep="\t", index=False, header=False)
    print("\n".join(snp_seq_list))


if __name__ == "__main__":
    typer.run(get_snp_sequence)
