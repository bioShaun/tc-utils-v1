from pathlib import Path

import pandas as pd
import typer
from loguru import logger


def make_vcf_header(chr_size: Path) -> str:
    header = "##fileformat=VCFv4.2\n"
    with chr_size.open() as f:
        for line in f:
            chr, size = line.strip().split("\t")[:2]
            header += f"##contig=<ID={chr},length={size}>\n"
        header += '#FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    return header


def main(gt_file: Path, chr_size: Path, vcf_file: Path):
    logger.info(f"loading gt file...")
    gt_df = pd.read_table(gt_file)
    gt_df["ID"] = gt_df["CHROM"] + "_" + gt_df["POS"].astype(str)
    gt_df["FILTER"] = "."
    gt_df["QUAL"] = "."
    gt_df["INFO"] = "."
    gt_df["FORMAT"] = "GT"
    vcf_df = gt_df[
        [
            "CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            *gt_df.columns[4:].tolist(),
        ]
    ].copy()
    vcf_df.rename(columns={"CHROM": "#CHROM"}, inplace=True)
    logger.info("make header...")
    header = make_vcf_header(chr_size)
    logger.info("make vcf...")
    with vcf_file.open("w") as f:
        f.write(header)
    vcf_df.to_csv(vcf_file, index=False, sep="\t", mode="a")


if __name__ == "__main__":
    typer.run(main)
