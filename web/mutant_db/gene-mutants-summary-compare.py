from pathlib import Path
from typing import List, Optional

import pandas as pd
import typer

HEADER_MAP = {
    "gene": "基因ID",
    "total_variant_count": "总变异数",
    "sample_count": "样本数",
    "sample_list": "样本列表",
    "shared_variant_count": "共有变异数",
    "stop_lost_gain_sample_list": "Stop Lost/Gain样本列表",
    "frame_shift_sample_list": "Frameshift样本列表",
    "missense_sample_list": "Missense样本列表",
}


def gene_mutants_summary(
    va_df: pd.DataFrame, sample_df: List[str], prefix: str
) -> pd.DataFrame:
    df = va_df[va_df["sample_id"].isin(sample_df)].copy()

    stop_df = df[(df["type"].str.contains("stop")) & (df["impact"] == "HIGH")]
    frame_shift_df = df[(df["type"].str.contains("frameshift"))]
    misssense_df = df[(df["type"].str.contains("missense"))]
    core_df = pd.concat([stop_df, frame_shift_df, misssense_df]).drop_duplicates()

    gene_sample_df = core_df.groupby(["gene"])["sample_id"].unique().reset_index()
    gene_sample_df["sample_count"] = gene_sample_df["sample_id"].map(len)
    gene_sample_df["sample_list"] = gene_sample_df["sample_id"].str.join(",")

    core_va_df = core_df.drop_duplicates(subset=["chrom", "pos", "refer", "alt"])
    va_count_df = (
        core_df.groupby(["chrom", "pos", "refer", "alt"])["sample_id"]
        .unique()
        .map(len)
        .reset_index()
    )
    shared_va_count_df = va_count_df[va_count_df["sample_id"] > 1]
    shared_core_df = core_va_df.merge(
        shared_va_count_df.drop(["sample_id"], axis=1),
        on=["chrom", "pos", "refer", "alt"],
        how="inner",
    )

    total_va_by_gene = (
        core_va_df["gene"].value_counts().rename("total_variant_count").reset_index()
    )
    shared_va_by_gene = (
        shared_core_df["gene"]
        .value_counts()
        .rename("shared_variant_count")
        .reset_index()
    )

    stop_gene_sample_df = stop_df.groupby(["gene"])["sample_id"].unique().reset_index()
    stop_gene_sample_df["stop_lost_gain_sample_list"] = stop_gene_sample_df[
        "sample_id"
    ].str.join(",")

    frame_shift_gene_sample_df = (
        frame_shift_df.groupby(["gene"])["sample_id"].unique().reset_index()
    )
    frame_shift_gene_sample_df["frame_shift_sample_list"] = frame_shift_gene_sample_df[
        "sample_id"
    ].str.join(",")

    misssense_gene_sample_df = (
        misssense_df.groupby(["gene"])["sample_id"].unique().reset_index()
    )
    misssense_gene_sample_df["missense_sample_list"] = misssense_gene_sample_df[
        "sample_id"
    ].str.join(",")

    merged_count = (
        total_va_by_gene.merge(shared_va_by_gene, how="left")
        .merge(gene_sample_df[["gene", "sample_count", "sample_list"]], how="left")
        .fillna(0)
    )
    merged_count_sample = (
        merged_count.merge(
            stop_gene_sample_df[["gene", "stop_lost_gain_sample_list"]], how="left"
        )
        .merge(
            frame_shift_gene_sample_df[["gene", "frame_shift_sample_list"]], how="left"
        )
        .merge(misssense_gene_sample_df[["gene", "missense_sample_list"]], how="left")
        .fillna("-")
    ).rename(columns=HEADER_MAP)
    merged_count_sample.columns = [
        "基因ID",
        *[f"{prefix}_{x}" for x in merged_count_sample.columns[1:]],
    ]
    return merged_count_sample


def main(
    db_file: Path, output_file: Path, target_sample_path: Path, all_sample_path: Path
):
    """Summarize gene mutations in a database file.

    Args:
        db_file (Path): Path to the database file.
        output_file (Path): Path to the output file.
        sample_list (Optional[Path]): Path to the sample list file.
    """
    # Read the database file
    df = pd.read_csv(db_file, sep="\t")

    # Filter the dataframe based on the sample list
    target_sample_list = pd.read_csv(target_sample_path, sep="\t", header=None)[
        0
    ].to_list()
    all_sample_list = pd.read_csv(all_sample_path, sep="\t", header=None)[0].to_list()
    background_sample_list = list(set(all_sample_list) - set(target_sample_list))
    # Summarize the gene mutations
    merged_count_sample = gene_mutants_summary(
        df, target_sample_list, "目标样品"
    ).merge(
        gene_mutants_summary(df, background_sample_list, "其他样品"),
        how="left",
    )

    int_cols = [each for each in merged_count_sample.columns if "数" in each]
    merged_count_sample[int_cols] = merged_count_sample[int_cols].fillna(0).astype(int)

    str_cols = [each for each in merged_count_sample.columns if "列表" in each]
    merged_count_sample[str_cols] = merged_count_sample[str_cols].fillna("-")

    merged_count_sample.to_excel(output_file, index=False)


if __name__ == "__main__":
    typer.run(main)
