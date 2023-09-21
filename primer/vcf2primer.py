import typer

from Bio import SeqIO
from pathlib import Path
import pandas as pd
from loguru import logger
from tqdm import tqdm


def main(vcf_file: Path, ref: Path, out_file: Path, flank_size: int = 200) -> None:
    ref_dict = {}
    logger.info("Loading reference")
    for record in SeqIO.parse(ref, format="fasta"):
        ref_dict[record.id] = record
    vcf_df = pd.read_table(
        vcf_file,
        comment="#",
        usecols=[0, 1, 3, 4],
        names=["chrom", "pos", "ref", "alt"],
    )
    out_list = []
    for row in tqdm(vcf_df.itertuples()):
        pos = row.pos - 1
        start = pos - 200 if pos >= 200 else 0
        left_seq = str(ref_dict[row.chrom][start : pos - 1])
        right_seq = str(ref_dict[row.chrom][pos + 1 : pos + 1 + 200])
        target_seq = str(ref_dict[row.chrom][pos])
        primer_seq = f"{left_seq}[{target_seq}/{row.alt}]{right_seq}"
        out_list.append({"name": f"{row.chrom}_{row.pos}", "sequence": primer_seq})
    out_df = pd.DataFrame(out_list)
    out_df.to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)
