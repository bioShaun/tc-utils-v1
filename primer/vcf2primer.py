from pathlib import Path

import pandas as pd
import typer
from Bio import SeqIO
from loguru import logger
from tqdm import tqdm


def main(vcf_file: Path, ref: Path, out_file: Path, flank_size: int = 200) -> None:
    vcf_df = pd.read_table(
        vcf_file,
        comment="#",
        usecols=[0, 1, 3, 4],
        names=["chrom", "pos", "ref", "alt"],
    )
    out_list = []

    for record in SeqIO.parse(ref, format="fasta"):
        logger.info(f"Process chrom: {record.id}")
        # ref_dict[record.id] = record
        vcf_df["chrom"] = vcf_df["chrom"].astype("str")
        chrom_vcf_df = vcf_df[vcf_df["chrom"] == record.id]
        for row in chrom_vcf_df.itertuples():
            pos = row.pos - 1
            start = pos - flank_size if pos >= flank_size else 0
            left_seq = str(record.seq[start:pos])
            right_seq = str(
                record.seq[pos + len(row.ref) : pos + len(row.ref) + flank_size]
            )
            # target_seq = str(record.seq[pos])
            primer_seq = f"{left_seq}[{row.ref}/{row.alt}]{right_seq}"
            out_list.append({"name": f"{row.chrom}_{row.pos}", "sequence": primer_seq})
    out_df = pd.DataFrame(out_list)
    out_df.to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)
