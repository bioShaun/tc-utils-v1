from pathlib import Path

import pandas as pd
import typer


def map_gt(gt: str) -> str:
    """Map a genotype string to its corresponding state."""
    if gt == "./.":
        return "MISS"
    ref, alt = gt.split("/")
    if ref == alt:
        return "HOM"
    return "HET"


def compare_gt(df: pd.DataFrame, compare_a: str, compare_b: str):
    stats_gt_df = df[[compare_a, compare_b]].map(map_gt)
    mask1 = stats_gt_df[compare_a] != "MISS"
    mask2 = stats_gt_df[compare_b] != "MISS"
    non_miss_df = stats_gt_df[mask1 & mask2]
    non_miss_count = len(non_miss_df)
    a_hom_count = len(non_miss_df[non_miss_df[compare_a] == "HOM"])
    a_het_count = len(non_miss_df[non_miss_df[compare_a] == "HET"])
    b_hom_count = len(non_miss_df[non_miss_df[compare_b] == "HOM"])
    b_het_count = len(non_miss_df[non_miss_df[compare_b] == "HET"])

    mask3 = non_miss_df[compare_a] == "HOM"
    mask4 = non_miss_df[compare_b] == "HOM"
    homo_non_miss_df = non_miss_df[mask3 & mask4]
    het_non_miss_df = non_miss_df[~non_miss_df.index.isin(homo_non_miss_df.index)]
    total_homo_count = len(homo_non_miss_df)
    total_het_count = len(het_non_miss_df)
    homo_compare_df = df[df.index.isin(homo_non_miss_df.index)]
    het_compare_df = df[df.index.isin(het_non_miss_df.index)]
    homo_equal_count = len(
        homo_compare_df[homo_compare_df[compare_a] == homo_compare_df[compare_b]]
    )
    homo_equal_count_out = (
        f"{100*homo_equal_count / total_homo_count:.2f}% ({homo_equal_count})"
    )
    het_equal_count = len(
        het_compare_df[het_compare_df[compare_a] == het_compare_df[compare_b]]
    )
    het_equal_count_out = (
        f"{100*het_equal_count / total_het_count:.2f}% ({het_equal_count})"
    )
    total_a_equal_b = homo_equal_count + het_equal_count
    total_a_equal_b_out = (
        f"{100*total_a_equal_b / non_miss_count:.2f}% ({total_a_equal_b})"
    )
    return {
        "A": compare_a,
        "B": compare_b,
        "总位点数": len(df),
        "有效位点": non_miss_count,
        "A_纯合": a_hom_count,
        "A_杂合": a_het_count,
        "B_纯合": b_hom_count,
        "B_杂合": b_het_count,
        "整体相似度": total_a_equal_b_out,
        "纯合相似度": homo_equal_count_out,
        "杂合相似度": het_equal_count_out,
    }


def main(genotype_file: Path, compare_list: Path, output_prefix: Path):
    df = pd.read_table(genotype_file)
    compare_df = pd.read_table(compare_list, header=None, names=["A", "B"])
    # gt_name_df = df.map(map_gt)
    # melt_gt_name_df = gt_name_df.melt(var_name="sample_id", value_name="genotype")
    # stats_count = (
    #     melt_gt_name_df.groupby("sample_id")["genotype"].value_counts().unstack(1)
    # )
    # stats_count["TOTAL"] = stats_count.sum(1)
    # stats_count["NON-MISS"] = stats_count["TOTAL"] - stats_count["MISS"]
    # stats_count.columns.name = ""
    # out_count_colums = ["NON-MISS", "HOM", "HET"]
    # out_count_colums_a = [f"A-{each}" for each in out_count_colums]
    # out_count_colums_b = [f"B-{each}" for each in out_count_colums]
    # stats_count_a = stats_count.reset_index().rename(columns={"sample_id": "A"})
    # stats_count_a.columns = out_count_colums_a
    # stats_count_b = stats_count.reset_index().rename(columns={"sample_id": "B"})
    # stats_count_b.columns = out_count_colums_b
    # compare_data_basic = compare_df.merge(stats_count_a).merge(stats_count_b)

    sample_stats_list = []
    chrom_stats_list = []
    for row in compare_df.itertuples():
        if row.A in df.columns and row.B in df.columns:
            sample_stats_list.append(compare_gt(df, row.A, row.B))
            for chrom, chrom_df in df.groupby("CHROM"):
                chrom_stats_dict = {"chrom": chrom}
                chrom_stats_dict.update(compare_gt(chrom_df, row.A, row.B))
                print(chrom_stats_dict)
                chrom_stats_list.append(chrom_stats_dict)
    sample_stats_df = pd.DataFrame(sample_stats_list)
    sample_stats = f"{output_prefix}.样品.xlsx"
    sample_stats_df.to_excel(sample_stats, index=False)
    sample_chrom_stats_df = pd.DataFrame(chrom_stats_list)
    sample_chrom_stats = f"{output_prefix}.染色体.xlsx"
    sample_chrom_stats_df.to_excel(sample_chrom_stats, index=False)


if __name__ == "__main__":
    typer.run(main)
