from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd
import typer
from tqdm import tqdm


@dataclass
class SiteStat:
    total_sites: int
    equal_sites: int
    ne_hom_sites: int
    ne_het_sites: int
    na_sites: int
    bc_rate: Optional[float]
    na_rate: Optional[float]
    het_rate: Optional[float]


def cal_bc_rate(equal_sites: int, ne_hom_sites: int, ne_het_sites: int) -> float:
    return (equal_sites + ne_het_sites * 0.5) / (
        equal_sites + ne_hom_sites + ne_het_sites
    )


def cal_reversion_rate(df: pd.DataFrame, parent_col: str, sample_col: str) -> SiteStat:
    if parent_col in df.columns:
        compare_df = df[[parent_col, sample_col]].copy()
        parent_hom_filter = compare_df[parent_col].isin(["0/0", "1/1"])
        total_sites = len(compare_df[parent_hom_filter])
        sample_non_missing_filter = compare_df[sample_col] != "./."
        compare_df = compare_df[parent_hom_filter & sample_non_missing_filter]
        equal_sites_count = len(
            compare_df[compare_df[sample_col] == compare_df[parent_col]]
        )
        ne_sites = compare_df[compare_df[sample_col] != compare_df[parent_col]]
        ne_hom_sites_count = len(ne_sites[ne_sites[sample_col].isin(["0/0", "1/1"])])
        ne_het_sites_count = len(ne_sites[ne_sites[sample_col].isin(["0/1", "1/0"])])
        bc_rate = (
            0
            if len(compare_df) == 0
            else cal_bc_rate(equal_sites_count, ne_hom_sites_count, ne_het_sites_count)
        )
    else:
        total_sites = len(df)
        compare_df = df[[sample_col]].copy()
        sample_non_missing_filter = compare_df[sample_col] != "./."
        compare_df = compare_df[sample_non_missing_filter]
        bc_rate = None
        equal_sites_count = 0
        ne_hom_sites_count = 0
    valid_sites = len(compare_df)
    na_sites = total_sites - valid_sites
    na_rate = None if total_sites == 0 else na_sites / total_sites
    het_count = len(compare_df[compare_df[sample_col] == "0/1"])
    het_rate = None if valid_sites == 0 else het_count / valid_sites

    return SiteStat(
        total_sites=total_sites,
        equal_sites=equal_sites_count,
        ne_hom_sites=ne_hom_sites_count,
        ne_het_sites=het_count,
        na_sites=na_sites,
        bc_rate=bc_rate,
        na_rate=na_rate,
        het_rate=het_rate,
    )


def main(
    gt_file: Path, parent_file: Path, output_file: Path, va_id: Optional[Path] = None
):
    gt = pd.read_csv(gt_file, sep="\t")
    if va_id is not None:
        va_id_list = pd.read_csv(va_id, sep="\t", header=None)[0].tolist()
        gt["variant_id"] = gt["CHROM"].astype(str) + "_" + gt["POS"].astype(str)
        gt = gt[gt["variant_id"].isin(va_id_list)]
        gt.drop("variant_id", axis=1, inplace=True)
    parent = pd.read_csv(parent_file, sep="\t")
    for row in tqdm(parent.itertuples(), total=len(parent)):
        parent_id = row.Parent_id
        sample_id = row.Sample_id
        site_stats = cal_reversion_rate(
            df=gt, parent_col=str(parent_id), sample_col=str(sample_id)
        )
        parent.loc[row.Index, "Total"] = site_stats.total_sites
        parent.loc[row.Index, "A"] = site_stats.equal_sites
        parent.loc[row.Index, "B"] = site_stats.ne_hom_sites
        parent.loc[row.Index, "H"] = site_stats.ne_het_sites
        parent.loc[row.Index, "N"] = site_stats.na_sites
        parent.loc[row.Index, "BC_rate"] = site_stats.bc_rate
        parent.loc[row.Index, "N_rate"] = site_stats.na_rate
        parent.loc[row.Index, "HET_rate"] = site_stats.het_rate
    parent = parent[parent["N_rate"] < 1]
    parent.to_excel(output_file, index=False, float_format="%.3f", na_rep="NA")


if __name__ == "__main__":
    typer.run(main)
