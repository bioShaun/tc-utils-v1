import gzip
from pathlib import Path

import delegator
import pandas as pd
import typer
from loguru import logger
from tqdm import tqdm


def main(
    vcf_file: Path, filter_vcf_file: Path, window_size: int = 120, snp_count: int = 4
) -> None:
    filter_vcf_info = open(filter_vcf_file, "w")
    with gzip.open(vcf_file, "rt") as f:
        for eachline in f:
            if eachline.startswith("#"):
                filter_vcf_info.write(eachline)
            else:
                break
    filter_vcf_info.close()

    dfs = pd.read_table(vcf_file, comment="#", header=None, chunksize=100_000)
    passed_dfs = []
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
        filter_df = df[~df.index.isin(cluster_row_indexes)]
        # passed_dfs.append(filter_df)
        logger.info(f"Writing {start_pos} - {end_pos}")
        filter_df.to_csv(filter_vcf_file, sep="\t", index=False, header=False, mode="a")
    logger.info(f"Compressing {filter_vcf_file}")
    delegator.run(
        f"/public/software/miniconda3/envs/plant-gmap/bin/bgzip -f {filter_vcf_file}"
    )


if __name__ == "__main__":
    typer.run(main)
