from pathlib import Path
from typing import List

import delegator
import pandas as pd
import typer
from loguru import logger


def vcf2depth(vcf_path: Path, force: bool = False) -> Path:
    dp_path = vcf_path.with_suffix(".dp.txt.gz")
    if dp_path.exists() and not force:
        return dp_path
    cmd = f"bcftools query -f '%CHROM\\t%POS\\t[\\t%DP]\n' ${vcf_path} | gzip > ${dp_path}"
    logger.info(cmd)
    delegator.run(cmd)
    return dp_path


def get_stats_df(df: pd.DataFrame) -> pd.DataFrame:
    min_cov = df.min(1)
    max_cov = df.max(1)
    ave_cov = df.sum(1) / df.shape[1]
    quantile_cov = df.quantile(0.5, axis=1)
    merged_df = pd.concat([min_cov, max_cov, ave_cov, quantile_cov], axis=1)
    merged_df.columns = ["min_cov", "max_cov", "mean_cov", "quantile_cov"]
    return merged_df


def main(
    vcf_path: Path,
    out_file: Path,
    cov: List[int] = [1, 5, 10, 20, 30, 50, 100],
):
    dp_path = vcf2depth(vcf_path)
    dp_df = pd.read_table(dp_path, header=None, index_col=[0, 1])
    dp_df = dp_df.applymap(lambda x: int(x) if str(x).isdigit() else 0)  # type: ignore
    stats_df = get_stats_df(dp_df)
    cov_df_list = []
    for cov_i in cov:
        logger.info(f"Calculate coverage >= {cov_i}x statistics ...")
        cov_i_df_matrix = dp_df >= cov_i
        cover_df = cov_i_df_matrix.sum(1)
        cover_ratio_df = cover_df / cov_i_df_matrix.shape[1]
        cover_ratio_df.name = f"coverage_{cov_i}x"
        cov_df_list.append(cover_ratio_df)
    cover_ratio_df = pd.concat([stats_df, *cov_df_list], axis=1)
    cover_ratio_df.index.names = ["chrom", "pos"]
    cover_ratio_df = cover_ratio_df.reset_index()
    if out_file.with_suffix(".tsv"):
        cover_ratio_df.to_csv(out_file, sep="\t", index=False)
    elif out_file.with_suffix(".xlsx"):
        cover_ratio_df.to_excel(out_file, index=False)
    else:
        logger.error(f"tsv or xlsx suffix required, {out_file.suffix} is not supported")


if __name__ == "__main__":
    typer.run(main)
