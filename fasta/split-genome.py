from pathlib import Path

import pandas as pd
import typer
from loguru import logger
from tqdm import tqdm

SPLIT_SIZE = 500_000_000

GFF_COLUMNS = [
    "seqid",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
]


def parse_gff_with_pandas(gff_file, feature_type=None):
    """使用pandas解析GFF文件"""
    # 读取GFF文件，跳过注释行
    df = pd.read_csv(
        gff_file,
        sep="\t",
        comment="#",
        names=GFF_COLUMNS,
        na_values=".",
        skip_blank_lines=True,
    )

    # 如果指定了feature类型，则过滤
    if feature_type:
        genes_df = df[df["feature"] == feature_type].copy()
    else:
        # 自动选择最合适的feature类型
        feature_counts = df["feature"].value_counts()

        # 优先级顺序
        priority_features = ["gene", "mRNA", "transcript", "CDS", "exon"]
        selected_feature = None

        for feat in priority_features:
            if feat in feature_counts.index and feature_counts[feat] > 0:
                selected_feature = feat
                break

        # 如果没有找到优先feature，使用数量最多的
        if selected_feature is None:
            selected_feature = feature_counts.index[0]

        genes_df = df[df["feature"] == selected_feature].copy()
        print(f"自动选择feature类型: {selected_feature} (共 {len(genes_df)} 个)")
        print(f"可用的feature类型: {', '.join(feature_counts.index[:10].tolist())}")

    # 确保坐标是整数类型
    genes_df["start"] = genes_df["start"].astype(int)
    genes_df["end"] = genes_df["end"].astype(int)

    return genes_df


def calculate_gaps(genes_df):
    """计算基因间隔"""
    # 按染色体和起始位置排序
    genes_sorted = genes_df.sort_values(["seqid", "start"]).reset_index(drop=True)

    # 创建间隔数据
    gaps = []

    # 按染色体分组
    for seqid, group in genes_sorted.groupby("seqid"):
        group = group.reset_index(drop=True)

        # 计算相邻基因间隔
        for i in tqdm(range(len(group) - 1), desc=f"计算染色体{seqid}的基因间隔"):
            gene1 = group.iloc[i]
            gene2 = group.iloc[i + 1]

            gap_size = gene2["start"] - gene1["end"] - 1

            if gap_size > 0:
                gaps.append(
                    {
                        "Chromosome": seqid,
                        "gap_size": gap_size,
                        "gap_start": gene1["end"] + 1,
                        "gap_end": gene2["start"] - 1,
                        "upstream_gene_end": gene1["end"],
                        "downstream_gene_start": gene2["start"],
                        "upstream_gene_strand": gene1["strand"],
                        "downstream_gene_strand": gene2["strand"],
                        "upstream_gene_info": gene1["attributes"],
                        "downstream_gene_info": gene2["attributes"],
                    }
                )

    # 转换为DataFrame
    gaps_df = pd.DataFrame(gaps)

    return gaps_df


def generate_split_chr_bed(best_split_site_df: pd.DataFrame) -> pd.DataFrame:
    split_df_list = []
    for row in best_split_site_df.itertuples():
        chrom = row.Chromosome
        start = row.gap_start
        end = row.gap_end
        split_point = (start + end) // 2
        split_df_list.append(
            {
                "Chromosome": chrom,
                "start": 0,
                "end": split_point,
                "id": f"{chrom}a",
            }
        )
        split_df_list.append(
            {
                "Chromosome": chrom,
                "start": split_point,
                "end": row.chrom_size,
                "id": f"{chrom}b",
            }
        )
    return pd.DataFrame(split_df_list)


def split_chrom(
    chr_df: pd.DataFrame, gap_df: pd.DataFrame, min_gene_gap: int
) -> pd.DataFrame:
    add_gap_df = gap_df.merge(chr_df)
    split_site_add_gap_df = add_gap_df[
        (add_gap_df["gap_start"] >= add_gap_df["split_start"])
        & (add_gap_df["gap_end"] <= add_gap_df["split_end"])
    ]
    best_split_site_idx = split_site_add_gap_df.groupby("Chromosome")[
        "gap_size"
    ].idxmax()

    best_split_site_df = split_site_add_gap_df.loc[best_split_site_idx].copy()
    # check if there is any chromosome gap size < min_gene_gap
    gap_size_not_passed_df = best_split_site_df[
        best_split_site_df["gap_size"] < min_gene_gap
    ]
    if gap_size_not_passed_df.empty:
        print("All chromosomes have split sites.")
        print("gene gap size 最小为:", best_split_site_df["gap_size"].min())
        return generate_split_chr_bed(best_split_site_df)
    raise ValueError(
        f"There are {len(not_in_best_split_site)} chromosomes: {not_in_best_split_site['Chromosome'].to_list()} not in best_split_site_df."
    )


def generate_split_gff(split_chr_bed: pd.DataFrame, gff: Path, out_gff: Path) -> None:
    gff_df = pd.read_table(
        gff,
        sep="\t",
        header=None,
        names=GFF_COLUMNS,
        skip_blank_lines=True,
        comment="#",
    )
    rename_split_chr_bed = split_chr_bed.copy()
    rename_split_chr_bed.columns = ["seqid", "split_start", "split_end", "new_seqid"]
    add_split_chr_df = gff_df.merge(
        rename_split_chr_bed,
    )
    filter1 = add_split_chr_df["start"] >= add_split_chr_df["split_start"]
    filter2 = add_split_chr_df["start"] <= add_split_chr_df["split_end"]
    filter_add_split_chr_df = add_split_chr_df[filter1 & filter2].copy()
    filter_add_split_chr_df["new_start"] = (
        filter_add_split_chr_df["start"] - filter_add_split_chr_df["split_start"]
    )
    filter_add_split_chr_df["new_end"] = (
        filter_add_split_chr_df["end"] - filter_add_split_chr_df["split_start"]
    )
    filter_add_split_chr_df.sort_values(
        ["new_seqid", "new_start", "new_end"], inplace=True
    )
    filter_add_split_chr_df.to_csv(
        out_gff,
        sep="\t",
        index=False,
        header=False,
        columns=[
            "new_seqid",
            "source",
            "feature",
            "new_start",
            "new_end",
            "score",
            "strand",
            "phase",
            "attributes",
        ],
    )


def main(
    genome_fai: Path, gff: Path, split_cat_bed: Path, min_gene_gap: int = 10_000
) -> None:
    fai_df = pd.read_table(
        genome_fai,
        header=None,
        names=["Chromosome", "chrom_size"],
        usecols=[0, 1],
        sep="\t",
    )
    fai_df["split_start"] = fai_df.apply(
        lambda x: max(0.3, 1 - SPLIT_SIZE / x["chrom_size"]) * x["chrom_size"], axis=1
    )
    fai_df["split_end"] = fai_df["chrom_size"] - fai_df["split_start"]
    need_to_split_df = fai_df[fai_df["chrom_size"] >= SPLIT_SIZE]
    do_not_need_to_split_df = fai_df[fai_df["chrom_size"] < SPLIT_SIZE].copy()

    gene_df = parse_gff_with_pandas(gff)
    gap_df = calculate_gaps(gene_df)
    split_chr_bed = split_chrom(need_to_split_df, gap_df, min_gene_gap)
    do_not_need_to_split_df["start"] = 0
    do_not_need_to_split_df["end"] = do_not_need_to_split_df["chrom_size"]
    do_not_need_to_split_df["id"] = do_not_need_to_split_df["Chromosome"]
    split_chr_bed_all = pd.concat([split_chr_bed, do_not_need_to_split_df])
    split_chr_bed_all.to_csv(
        split_cat_bed,
        sep="\t",
        index=False,
        header=False,
        columns=split_chr_bed.columns,
    )
    out_gff = gff.with_suffix(".split.gff")
    logger.info(f"生成分割后的GFF文件: {out_gff}")
    generate_split_gff(split_chr_bed, gff, out_gff)


if __name__ == "__main__":
    typer.run(main)
