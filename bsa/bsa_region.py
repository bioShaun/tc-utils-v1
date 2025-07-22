#!/usr/bin/env python3
"""
BSA (Bulked Segregant Analysis) 阈值片段识别脚本
用于识别QTL-seq或BSA分析中超过阈值的连续片段
"""

import argparse
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


class BSARegionFinder:
    def __init__(
        self,
        threshold: float,
        min_snps: int = 5,
        min_length: int = 100000,
        max_gap: int = 500000,
        window_size: int = 1000000,
    ):
        """
        初始化BSA区间识别器

        Parameters:
        -----------
        threshold : float
            统计指标的阈值（如Δ(SNP-index)、G'值等）
        min_snps : int
            最小SNP数量
        min_length : int
            最小片段长度（bp）
        max_gap : int
            最大允许gap距离（bp）
        window_size : int
            滑动窗口大小（bp）
        """
        self.threshold = threshold
        self.min_snps = min_snps
        self.min_length = min_length
        self.max_gap = max_gap
        self.window_size = window_size

    def load_data(
        self,
        file_path: str,
        chr_col: str = "CHROM",
        pos_col: str = "POS",
        stat_col: str = "delta_snp_index",
    ) -> pd.DataFrame:
        """
        加载BSA分析数据

        Parameters:
        -----------
        file_path : str
            输入文件路径
        chr_col : str
            染色体列名
        pos_col : str
            位置列名
        stat_col : str
            统计指标列名
        """
        try:
            # 尝试读取不同格式的文件
            if file_path.endswith(".csv"):
                df = pd.read_csv(file_path)
            elif file_path.endswith((".txt", ".tsv")):
                df = pd.read_csv(file_path, sep="\t")
            else:
                # 尝试自动检测分隔符
                df = pd.read_csv(file_path, sep=None, engine="python")

            # 检查必需的列
            required_cols = [chr_col, pos_col, stat_col]
            missing_cols = [col for col in required_cols if col not in df.columns]
            if missing_cols:
                raise ValueError(f"缺少必需的列: {missing_cols}")

            # 数据预处理
            df = df.dropna(subset=[chr_col, pos_col, stat_col])
            df[pos_col] = df[pos_col].astype(int)
            df = df.sort_values([chr_col, pos_col])

            print(f"成功加载数据: {len(df)} 行，{len(df.columns)} 列")
            return df

        except Exception as e:
            print(f"加载数据失败: {e}")
            sys.exit(1)

    def find_threshold_regions(
        self,
        df: pd.DataFrame,
        chr_col: str = "CHROM",
        pos_col: str = "POS",
        stat_col: str = "delta_snp_index",
        use_abs: bool = True,
    ) -> List[Dict]:
        """
        识别超过阈值的连续片段

        Parameters:
        -----------
        df : pd.DataFrame
            输入数据
        chr_col : str
            染色体列名
        pos_col : str
            位置列名
        stat_col : str
            统计指标列名
        use_abs : bool
            是否使用绝对值判断阈值
        """
        regions = []

        for chrom in sorted(df[chr_col].unique()):
            chr_data = df[df[chr_col] == chrom].copy()
            chr_data = chr_data.sort_values(pos_col)

            # 判断是否超过阈值
            if use_abs:
                above_threshold = abs(chr_data[stat_col]) >= self.threshold
            else:
                above_threshold = chr_data[stat_col] >= self.threshold

            chr_data["above_threshold"] = above_threshold

            # 找到连续的超阈值区间
            chr_regions = self._find_continuous_regions(chr_data, chrom, pos_col)
            regions.extend(chr_regions)

        return regions

    def _find_continuous_regions(
        self, chr_data: pd.DataFrame, chrom: str, pos_col: str
    ) -> List[Dict]:
        """在单个染色体上找到连续区间"""
        regions = []
        current_region_start = None
        current_region_snps = []

        for idx, row in chr_data.iterrows():
            pos = row[pos_col]
            above_threshold = row["above_threshold"]

            if above_threshold:
                if current_region_start is None:
                    # 开始新区间
                    current_region_start = pos
                    current_region_snps = [idx]
                else:
                    # 检查是否在允许的gap范围内
                    last_pos = chr_data.loc[current_region_snps[-1], pos_col]
                    if pos - last_pos <= self.max_gap:
                        current_region_snps.append(idx)
                    else:
                        # gap太大，结束当前区间，开始新区间
                        self._finalize_region(
                            chr_data,
                            chrom,
                            pos_col,
                            current_region_start,
                            current_region_snps,
                            regions,
                        )
                        current_region_start = pos
                        current_region_snps = [idx]
            else:
                if current_region_start is not None:
                    # 结束当前区间
                    self._finalize_region(
                        chr_data,
                        chrom,
                        pos_col,
                        current_region_start,
                        current_region_snps,
                        regions,
                    )
                    current_region_start = None
                    current_region_snps = []

        # 处理最后一个区间
        if current_region_start is not None:
            self._finalize_region(
                chr_data,
                chrom,
                pos_col,
                current_region_start,
                current_region_snps,
                regions,
            )

        return regions

    def _finalize_region(
        self,
        chr_data: pd.DataFrame,
        chrom: str,
        pos_col: str,
        start_pos: int,
        snp_indices: List,
        regions: List[Dict],
    ):
        """完成区间并添加到结果中"""
        if len(snp_indices) < self.min_snps:
            return

        end_pos = chr_data.loc[snp_indices[-1], pos_col]
        length = end_pos - start_pos + 1

        if length < self.min_length:
            return

        region = {
            "chr": chrom,
            "start": start_pos,
            "end": end_pos,
            "length": length,
            "no_snps": len(snp_indices),
        }

        regions.append(region)

    def annotate_genes(
        self, regions: List[Dict], gene_file: Optional[str] = None
    ) -> List[Dict]:
        """
        为区间添加基因注释信息

        Parameters:
        -----------
        regions : List[Dict]
            识别的区间列表
        gene_file : str, optional
            基因注释文件路径 (GFF/GTF格式)
        """
        if gene_file is None:
            # 如果没有基因文件，设置no.genes为NA
            for region in regions:
                region["no.genes"] = "NA"
            return regions

        try:
            # 简单的基因计数（这里需要根据实际的基因注释格式调整）
            gene_df = pd.read_csv(
                gene_file,
                sep="\t",
                comment="#",
                header=None,
                names=[
                    "chr",
                    "source",
                    "feature",
                    "start",
                    "end",
                    "score",
                    "strand",
                    "frame",
                    "attributes",
                ],
            )

            # 只保留基因记录
            gene_df = gene_df[gene_df["feature"] == "gene"]

            for region in regions:
                # 计算重叠的基因数量
                chr_genes = gene_df[gene_df["chr"] == region["chr"]]
                overlapping_genes = chr_genes[
                    (chr_genes["start"] <= region["end"])
                    & (chr_genes["end"] >= region["start"])
                ]
                region["no.genes"] = len(overlapping_genes)

        except Exception as e:
            print(f"基因注释失败: {e}")
            for region in regions:
                region["no.genes"] = "NA"

        return regions

    def save_results(self, regions: List[Dict], output_file: str):
        """保存结果到文件"""
        if not regions:
            print("未找到符合条件的区间")
            return

        results_df = pd.DataFrame(regions)

        # 重新排列列的顺序
        column_order = ["chr", "start", "end", "length", "no.genes"]
        if "no_snps" in results_df.columns:
            column_order.append("no_snps")

        results_df = results_df[column_order]

        # 保存到文件
        results_df.to_csv(output_file, sep="\t", index=False)
        print(f"结果已保存到: {output_file}")
        print(f"共找到 {len(regions)} 个区间")

        # 打印前几行作为预览
        print("\n前5行结果预览:")
        print(results_df.head().to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description="BSA阈值片段识别")
    parser.add_argument("input_file", help="输入数据文件")
    parser.add_argument("-o", "--output", default="bsa_regions.txt", help="输出文件名")
    parser.add_argument(
        "-t", "--threshold", type=float, default=0.5, help="阈值 (默认: 0.5)"
    )
    parser.add_argument("--chr-col", default="CHROM", help="染色体列名")
    parser.add_argument("--pos-col", default="POS", help="位置列名")
    parser.add_argument("--stat-col", default="delta_snp_index", help="统计指标列名")
    parser.add_argument("--min-snps", type=int, default=5, help="最小SNP数量")
    parser.add_argument(
        "--min-length", type=int, default=100000, help="最小片段长度(bp)"
    )
    parser.add_argument("--max-gap", type=int, default=500000, help="最大gap距离(bp)")
    parser.add_argument("--gene-file", help="基因注释文件(GFF/GTF)")
    parser.add_argument("--no-abs", action="store_true", help="不使用绝对值判断阈值")

    args = parser.parse_args()

    # 创建识别器
    finder = BSARegionFinder(
        threshold=args.threshold,
        min_snps=args.min_snps,
        min_length=args.min_length,
        max_gap=args.max_gap,
    )

    # 加载数据
    df = finder.load_data(args.input_file, args.chr_col, args.pos_col, args.stat_col)

    # 找到阈值区间
    regions = finder.find_threshold_regions(
        df, args.chr_col, args.pos_col, args.stat_col, use_abs=not args.no_abs
    )

    # 基因注释
    regions = finder.annotate_genes(regions, args.gene_file)

    # 保存结果
    finder.save_results(regions, args.output)


if __name__ == "__main__":
    # 示例用法
    if len(sys.argv) == 1:
        print("BSA阈值片段识别脚本")
        print("\n基本用法:")
        print("python bsa_regions.py input_data.txt")
        print("\n完整用法:")
        print(
            "python bsa_regions.py input_data.txt -o output.txt -t 0.6 --min-snps 10 --gene-file genes.gff"
        )
        print("\n参数说明:")
        print("  input_file: 输入数据文件(CSV/TSV格式)")
        print("  -t: 阈值(默认0.5)")
        print("  --min-snps: 最小SNP数量(默认5)")
        print("  --min-length: 最小片段长度bp(默认100000)")
        print("  --max-gap: 最大gap距离bp(默认500000)")
        print("  --gene-file: 基因注释文件")
    else:
        main()
