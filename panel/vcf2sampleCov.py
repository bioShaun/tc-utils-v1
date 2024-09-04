from enum import Enum
from pathlib import Path
from typing import List, Tuple

import delegator
import pandas as pd
import typer
from loguru import logger


class OutType(str, Enum):
    vcf = "vcf"
    xlsx = "xlsx"


def vcf2depth(vcf_path: Path, force: bool = False) -> Path:
    dp_path = vcf_path.with_suffix(".dp.txt.gz")
    if dp_path.exists() and not force:
        return dp_path
    cmd = f"bcftools query -f '%CHROM\\t%POS[\\t%DP]\\n' {vcf_path} | gzip > {dp_path}"
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


def load_dp_df(vcf_path: Path, sample_file: Path) -> pd.DataFrame:
    dp_path = vcf2depth(vcf_path)
    sample_list = pd.read_csv(sample_file, header=None)[0].tolist()
    dp_df = pd.read_table(
        dp_path, header=None, index_col=[0, 1], names=["CHROM", "POS", *sample_list]
    )
    return dp_df.applymap(lambda x: int(x) if str(x).isdigit() else 0)  # type: ignore


def sample_site_coverage(
    dp_df: pd.DataFrame,
    cov: List[int] = [1, 5, 10, 20, 30, 50, 100],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    stats_df = get_stats_df(dp_df)
    site_cov_df_list = []
    sample_cov_df_list = []
    for cov_i in cov:
        logger.info(f"Calculate coverage >= {cov_i}x statistics ...")
        cov_i_df_matrix = dp_df >= cov_i
        site_cover_df = cov_i_df_matrix.sum(1)
        site_cover_ratio_df = site_cover_df / cov_i_df_matrix.shape[1]
        site_cover_ratio_df.name = f"coverage_{cov_i}x"
        site_cov_df_list.append(site_cover_ratio_df)
        sample_cover_df = cov_i_df_matrix.sum()
        sample_cover_ratio_df = sample_cover_df / cov_i_df_matrix.shape[0]
        sample_cover_ratio_df.name = f"coverage_{cov_i}x"
        sample_cov_df_list.append(sample_cover_ratio_df)
    sample_mean_cov_df = pd.DataFrame(dp_df.mean())
    sample_mean_cov_df.columns = ["mean_coverage"]
    cover_ratio_df = pd.concat([stats_df, *site_cov_df_list], axis=1)
    cover_ratio_df.index.names = ["chrom", "pos"]
    cover_ratio_df = cover_ratio_df.reset_index()
    sample_ratio_df = pd.concat([sample_mean_cov_df, *sample_cov_df_list], axis=1)
    sample_ratio_df.index.names = ["sample_id"]
    sample_ratio_df = sample_ratio_df.reset_index()
    return cover_ratio_df, sample_ratio_df


def main(
    vcf_path: Path,
    sample_file: Path,
    out_file_prefix: Path,
    out_type: OutType = OutType.xlsx,
    cov: List[int] = [1, 5, 10, 20, 30, 50, 100],
) -> None:
    dp_df = load_dp_df(vcf_path, sample_file)

    site_cover_ratio_df, sample_cover_ratio_df = sample_site_coverage(dp_df, cov)
    site_cover_ratio_file = f"{out_file_prefix}.site.coverage.{out_type.value}"
    sample_cover_ratio_file = f"{out_file_prefix}.sample.coverage.{out_type.value}"
    if out_type == "tsv":
        site_cover_ratio_df.to_csv(site_cover_ratio_file, sep="\t", index=False)
        sample_cover_ratio_df.to_csv(sample_cover_ratio_file, sep="\t", index=False)
    else:
        site_cover_ratio_df.to_excel(site_cover_ratio_file, index=False)
        sample_cover_ratio_df.to_excel(sample_cover_ratio_file, index=False)


if __name__ == "__main__":
    typer.run(main)
