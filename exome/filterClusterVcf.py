import gzip
from pathlib import Path

import delegator
import pandas as pd
import typer
from loguru import logger
from tqdm import tqdm


def filter_vcf(
    vcf_file: Path,
    output_file: Path,
    window_size: int = 60,
    snp_count: int = 2,
    bgzip_bin_path: str = "/public/home/zxchen/software/miniconda3/envs/exome/bin/",
) -> None:
    with open(output_file, "w") as output_file_handle:
        with gzip.open(vcf_file, "rt") as vcf_file_handle:
            for line in vcf_file_handle:
                if line.startswith("#"):
                    output_file_handle.write(line)
                else:
                    break
    dfs = pd.read_table(vcf_file, comment="#", header=None, chunksize=100_000)
    for df in dfs:
        cluster_rows = set()
        compare_rows = []
        start_pos = f"{df.iloc[0][0]}-{df.iloc[0][1]}"
        end_pos = f"{df.iloc[-1][0]}-{df.iloc[-1][1]}"
        logger.info(f"Processing {start_pos} - {end_pos}")
        for row in tqdm(df.itertuples()):
            if len(compare_rows) <= snp_count:
                compare_rows.append([row[0], row[1], row[2]])
                continue
            compare_row = compare_rows[1 - snp_count]
            if row[1] == compare_row[1] and row[2] - compare_row[2] <= window_size:
                for each_row in compare_rows:
                    cluster_rows.add(each_row[0])
                cluster_rows.add(row[0])
            compare_rows.append([row[0], row[1], row[2]])
            compare_rows = compare_rows[1:]
        cluster_row_indexes = list(cluster_rows)
        filtered_df = df[~df.index.isin(cluster_row_indexes)]
        logger.info(f"Writing {start_pos} - {end_pos}")
        filtered_df.to_csv(output_file, sep="\t", index=False, header=False, mode="a")
    logger.info(f"Compressing {output_file}")
    delegator.run(f"{bgzip_bin_path}/bgzip -f {output_file}")


if __name__ == "__main__":
    typer.run(filter_vcf)
