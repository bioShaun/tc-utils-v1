from pathlib import Path
import pandas as pd
import typer


BASE_DIR = "/data2/app/tcuni-project-system/server/public/analysis/BSA/"


def main(
    bsa_id: int, chrom: str, start: int, end: int, snp_file: Path, unit: str = "M"
) -> None:
    df_list = []
    bsa_dir = Path(BASE_DIR) / f"{bsa_id}/BSA-results"
    if unit == "M":
        start = start * 1_000_000
        end = end * 1_000_000
    for data_i in bsa_dir.glob(f"*/data/{chrom}.*"):
        if f"{chrom}.ed" in data_i.name:
            continue
        df = pd.read_csv(data_i)
        vcf_df = df[(df["POS"] >= start) & (df["POS"] <= end)].copy()
        print(vcf_df)
        vcf_df = vcf_df[vcf_df["Type"] == "SNP"].copy()
        vcf_df = vcf_df[["CHROM", "POS", "REF", "ALT"]].copy()
        vcf_df.drop_duplicates(inplace=True)
        vcf_df["ID"] = "."
        df_list.append(vcf_df[["CHROM", "POS", "ID", "REF", "ALT"]])
    merged_df = pd.concat(df_list)
    merged_df.to_csv(snp_file, sep="\t", header=False, index=False)


if __name__ == "__main__":
    typer.run(main)
