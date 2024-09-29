from io import StringIO
from pathlib import Path
from typing import Tuple

import pandas as pd
import typer


def flat_gt_table(gt_table: Path) -> pd.DataFrame:
    flat_data_list = []
    with open(gt_table) as gt_inf:
        for line in gt_inf:
            chrom, pos, refer, alt, *score_and_sample = line.strip().split("\t")
            try:
                for sample_i in score_and_sample[1:]:
                    flat_data_list.append(f"{chrom}\t{pos}\t{refer}\t{alt}\t{sample_i}")
            except ValueError:
                print(line)
                raise ValueError("sample list wrong type")
    flat_data_str = "\n".join(flat_data_list)
    return pd.read_table(
        StringIO(flat_data_str),
        header=None,
        names=["chrom", "pos", "refer", "alt", "accession"],
    )


def gt_from_af(af: float) -> str:
    if af < 0.25:
        return "0/0"
    elif 0.25 <= af < 0.75:
        return "0/1"
    return "1/1"


def split_accession(accession: str) -> Tuple[str, int, int, int, float, str]:
    sample_id, allele_info = accession.split("=")
    genotype, allele_depth = allele_info.split(";")
    if "." in genotype:
        genotype = "./."
        reference_depth, alternate_depth = 0, 0
    else:
        allele_depth = allele_depth.replace(".", "0")
        reference_depth, alternate_depth = [
            int(each) for each in allele_depth.split(",")
        ]
    total_depth = reference_depth + alternate_depth
    allele_freq = 0 if total_depth == 0 else alternate_depth / total_depth
    genotype = gt_from_af(allele_freq)
    return (
        sample_id,
        int(reference_depth),
        int(alternate_depth),
        total_depth,
        allele_freq,
        genotype,
    )


def ems_change(row: pd.Series) -> str:
    dna_change = f"{row.refer}-{row.alt}"
    if dna_change not in ["C-T", "G-A"]:
        dna_change = "other"
    return dna_change


def main(gt_table: Path) -> None:
    flat_gt_df = flat_gt_table(gt_table)
    flat_gt_df[
        [
            "sample_id",
            "ref_depth",
            "alt_depth",
            "total_depth",
            "allele_freq",
            "raw_genotype",
        ]
    ] = (
        flat_gt_df["accession"].map(split_accession).to_list()  # type: ignore
    )

    filter_flat_gt_df = flat_gt_df[flat_gt_df["raw_genotype"] != "0/0"].copy()
    filter_flat_gt_df["genotype"] = filter_flat_gt_df["raw_genotype"].map(
        lambda x: "HET" if x == "0/1" else "HOM"
    )
    filter_flat_gt_df["ems_change"] = filter_flat_gt_df.apply(ems_change, axis=1)

    stats_count = (
        filter_flat_gt_df.groupby(["sample_id", "genotype", "ems_change"])
        .size()
        .unstack(2)
        .fillna(0)
        .astype("int")
    )
    stats_ratio = (stats_count.T / stats_count.sum(1)).T
    stats_ratio = stats_ratio.map(lambda x: f"{x:.2%}")
    stats_ratio.columns = [f"{each}-Ratio" for each in stats_ratio.columns]
    stats_count.columns = [f"{each}-Count" for each in stats_count.columns]
    merged_stats = stats_count.merge(stats_ratio, left_index=True, right_index=True)
    merged_stats = merged_stats.reset_index()

    all_stats_count = (
        filter_flat_gt_df.groupby(["sample_id", "ems_change"])
        .size()
        .unstack(1)
        .fillna(0)
        .astype("int")
    )
    all_stats_ratio = (all_stats_count.T / all_stats_count.sum(1)).T
    all_stats_ratio = all_stats_ratio.map(lambda x: f"{x:.2%}")
    all_stats_ratio.columns = [f"{each}-Ratio" for each in all_stats_ratio.columns]
    all_stats_count.columns = [f"{each}-Count" for each in all_stats_count.columns]
    all_merged_stats = all_stats_count.merge(
        all_stats_ratio, left_index=True, right_index=True
    )  # type: ignore
    all_merged_stats = all_merged_stats.reset_index()
    all_merged_stats["genotype"] = "ALL"

    hom_het_all_merged_stats = pd.concat([merged_stats, all_merged_stats])
    hom_het_all_merged_stats.sort_values(["sample_id", "genotype"], inplace=True)

    stats_file = gt_table.with_suffix(".stats.xlsx")
    hom_het_all_merged_stats.to_csv(stats_file, index=False)


if __name__ == "__main__":
    typer.run(main)
