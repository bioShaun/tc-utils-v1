from pathlib import Path

import pandas as pd
import typer


def main(bsr_dir: Path, snp_file: Path, snp_only: bool = False) -> None:
    df_list = []
    for data_i in bsr_dir.glob(f"data/*.gz"):
        if "ed.refFreq" in data_i.name:
            continue
        vcf_df = pd.read_csv(data_i)
        if snp_only:
            vcf_df = vcf_df[vcf_df["Type"] == "SNP"].copy()
        vcf_df = vcf_df[["CHROM", "POS", "REF", "ALT"]].copy()
        vcf_df.drop_duplicates(inplace=True)
        vcf_df["ID"] = "."
        df_list.append(vcf_df[["CHROM", "POS", "ID", "REF", "ALT"]])
    merged_df = pd.concat(df_list)
    merged_df.to_csv(snp_file, sep="\t", header=False, index=False)


if __name__ == "__main__":
    typer.run(main)
