import gzip
from pathlib import Path

import delegator
import pandas as pd
import typer
from loguru import logger
from tqdm import tqdm


def main(
    csv_file: Path, filter_vcf_file: Path, window_size: int = 120, snp_count: int = 4
) -> None:
    df = pd.read_csv(
        csv_file, header=None, names=["chrom", "pos", "alleles", "missing", "het", "af"]
    )
    df.sort_values(by=["chrom", "pos"], inplace=True)
    cluster_rows = set()
    compare_rows = []
    for row in tqdm(df.itertuples()):
        if len(compare_rows) < snp_count - 1:
            compare_rows.append([row.Index, row.chrom, row.pos])
            continue
        compare_row = compare_rows[1 - snp_count]
        if row.chrom == compare_row[1] and row.pos - compare_row[2] <= window_size:
            for each_row in compare_rows:
                cluster_rows.add(each_row[0])
            cluster_rows.add(row.Index)
        compare_rows.append([row.Index, row.chrom, row.pos])
        compare_rows = compare_rows[1:]
    cluster_row_indexes = list(cluster_rows)
    filter_df = df[~df.index.isin(cluster_row_indexes)]
    # passed_dfs.append(filter_df)
    filter_df.to_csv(filter_vcf_file, index=False, header=False)


if __name__ == "__main__":
    typer.run(main)
