import gzip
import logging
import os

import allel
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import typer
from tqdm import tqdm

# 设置日志
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)
from pathlib import Path


def hdplot_large_vcf(vcf_file, chunk_size=1000, output_prefix=None, temp_dir=None):
    """
    为大型VCF文件优化的HDPlot实现，使用分块处理策略

    参数:
    vcf_file -- VCF文件路径
    chunk_size -- 每次处理的变异位点数量
    output_prefix -- 输出文件前缀
    temp_dir -- 临时文件目录

    返回:
    包含HDPlot统计量的DataFrame
    """
    # 创建临时目录
    if temp_dir is None:
        temp_dir = "hdplot_temp"
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # 检查VCF文件
    if not os.path.exists(vcf_file):
        raise FileNotFoundError(f"找不到VCF文件: {vcf_file}")

    logger.info(f"开始处理VCF文件: {vcf_file}")

    # 获取VCF文件的总变异位点数和样本数
    total_variants, samples = get_vcf_info(vcf_file)
    logger.info(f"VCF文件包含 {total_variants} 个变异位点和 {len(samples)} 个样本")

    # 计算需要处理的块数
    num_chunks = (total_variants + chunk_size - 1) // chunk_size
    logger.info(f"将分 {num_chunks} 块处理数据，每块 {chunk_size} 个变异位点")

    # 分块处理VCF文件
    chunk_results = []
    for chunk_idx in tqdm(range(num_chunks), desc="处理VCF块"):
        start_pos = chunk_idx * chunk_size
        # 处理当前块
        chunk_result = process_vcf_chunk(vcf_file, start_pos, chunk_size)

        # 保存中间结果
        if output_prefix:
            chunk_file = os.path.join(
                temp_dir, f"{output_prefix}_chunk_{chunk_idx}.csv"
            )
            chunk_result.to_csv(chunk_file, index=False)

        chunk_results.append(chunk_result)

        # 释放内存
        del chunk_result

    # 合并所有块的结果
    logger.info("合并所有块的结果")
    if output_prefix:
        # 从保存的文件中读取结果
        all_results = pd.concat(
            [
                pd.read_csv(os.path.join(temp_dir, f"{output_prefix}_chunk_{i}.csv"))
                for i in range(num_chunks)
            ]
        )
    else:
        # 直接合并内存中的结果
        all_results = pd.concat(chunk_results)

    # 保存最终结果
    if output_prefix:
        final_output = f"{output_prefix}_hdplot_results.csv"
        all_results.to_csv(final_output, index=False)
        logger.info(f"HDPlot结果已保存至: {final_output}")

        # 清理临时文件
        if os.path.exists(temp_dir):
            for f in os.listdir(temp_dir):
                if f.startswith(f"{output_prefix}_chunk_"):
                    os.remove(os.path.join(temp_dir, f))
            os.rmdir(temp_dir)

    return all_results


def get_vcf_info(vcf_file):
    """获取VCF文件的基本信息"""
    # 使用pysam或cyvcf2等库可以更高效地获取这些信息
    # 这里使用allel简单实现
    try:
        # 只读取头部信息
        vcf = allel.read_vcf_headers(vcf_file)
        samples = vcf.samples

        # 计算变异位点数量
        # 注意：这可能会消耗大量内存，实际应用中应使用更高效的方法
        # 如使用tabix索引或直接解析VCF文件
        if vcf_file.endswith(".gz"):
            with gzip.open(vcf_file, "rt") as f:
                total_variants = sum(1 for line in f if not line.startswith("#"))
        else:
            with open(vcf_file, "r") as f:
                total_variants = sum(1 for line in f if not line.startswith("#"))

        return total_variants, samples
    except Exception as e:
        logger.error(f"获取VCF信息失败: {str(e)}")
        raise


def process_vcf_chunk(vcf_file, start_pos, chunk_size):
    """处理VCF文件的一个块"""
    try:
        # 使用region参数或其他方法只加载特定范围的变异位点
        # 注意：allel.read_vcf不直接支持按位置范围读取
        # 这里使用一个简化的方法，实际应用中可能需要更复杂的实现

        # 使用tabix索引的方法（如果VCF已索引）或其他方法
        # 这里简化为读取整个文件然后切片，实际应用中不建议这样做
        callset = allel.read_vcf(
            vcf_file,
            region=None,
            fields=[
                "variants/CHROM",
                "variants/POS",
                "variants/ID",
                "calldata/GT",
                "calldata/AD",
            ],
        )

        # 切片获取当前块的数据
        end_pos = min(start_pos + chunk_size, len(callset["variants/POS"]))

        # 提取必要的信息
        chrom = callset["variants/CHROM"][start_pos:end_pos]
        pos = callset["variants/POS"][start_pos:end_pos]
        id_field = callset.get("variants/ID", np.repeat(".", len(pos)))[
            start_pos:end_pos
        ]
        genotypes = allel.GenotypeArray(callset["calldata/GT"][start_pos:end_pos])

        # 提取AD字段(等位基因深度)
        if "calldata/AD" in callset:
            ad_field = callset["calldata/AD"][start_pos:end_pos]
        else:
            raise ValueError("VCF文件中缺少AD字段，无法计算等位基因深度")

        # 计算HDPlot统计量
        return calculate_hdplot_stats(chrom, pos, id_field, genotypes, ad_field)

    except Exception as e:
        logger.error(f"处理VCF块失败: {str(e)}")
        raise


def calculate_hdplot_stats(chrom, pos, id_field, genotypes, ad_field):
    """计算HDPlot统计量"""
    n_variants = len(pos)
    n_samples = genotypes.shape[1]

    # 初始化结果DataFrame
    hdplot_table = pd.DataFrame(
        {
            "CHROM": chrom,
            "POS": pos,
            "ID": id_field,
            "depth_a": np.zeros(n_variants),
            "depth_b": np.zeros(n_variants),
            "ratio": np.zeros(n_variants),
            "num_hets": np.zeros(n_variants),
            "num_samples": np.repeat(n_samples, n_variants),
            "num_called": np.zeros(n_variants),
            "H_all": np.zeros(n_variants),
            "H": np.zeros(n_variants),
            "std": np.zeros(n_variants),
            "D": np.zeros(n_variants),
        }
    )

    # 处理基因型数据
    het_matrix = np.zeros((n_variants, n_samples), dtype=int)
    called_genos = np.zeros((n_variants, n_samples), dtype=int)

    # 使用numpy向量化操作提高性能
    is_het = genotypes.is_het()
    is_called = ~genotypes.is_missing()

    het_matrix = is_het.astype(int)
    called_genos = is_called.astype(int)

    # 提取等位基因读数
    allele_reads_1 = ad_field[:, :, 0]
    allele_reads_2 = ad_field[:, :, 1]

    # 只保留杂合子的读数
    allele_reads_1_het = allele_reads_1 * het_matrix
    allele_reads_2_het = allele_reads_2 * het_matrix

    # 计算每个位点的总读数和比例
    a_reads = np.sum(allele_reads_1_het, axis=1)
    b_reads = np.sum(allele_reads_2_het, axis=1)
    total_reads = a_reads + b_reads

    # 避免除以零
    ratio = np.zeros(n_variants)
    valid_idx = total_reads > 0
    ratio[valid_idx] = a_reads[valid_idx] / total_reads[valid_idx]

    # 计算标准差和Z分数
    std = np.sqrt(total_reads * 0.5 * 0.5)
    z = np.zeros(n_variants)
    valid_idx = std > 0
    z[valid_idx] = -((total_reads[valid_idx] / 2) - a_reads[valid_idx]) / std[valid_idx]

    # 计算杂合度
    num_hets = np.sum(het_matrix, axis=1)
    num_genos = np.sum(called_genos, axis=1)
    het_perc = num_hets / n_samples
    h = np.zeros(n_variants)
    valid_idx = num_genos > 0
    h[valid_idx] = num_hets[valid_idx] / num_genos[valid_idx]

    # 填充结果表
    hdplot_table["depth_a"] = a_reads
    hdplot_table["depth_b"] = b_reads
    hdplot_table["ratio"] = ratio
    hdplot_table["num_hets"] = num_hets
    hdplot_table["num_called"] = num_genos
    hdplot_table["H_all"] = het_perc
    hdplot_table["H"] = h
    hdplot_table["std"] = std
    hdplot_table["D"] = z

    return hdplot_table


def plot_hdplot(
    hdplot_table, output_file=None, threshold_h=0.5, threshold_d=3, max_points=10000
):
    """
    绘制HDPlot散点图，针对大型数据集进行了优化

    参数:
    hdplot_table -- HDPlot函数返回的DataFrame
    output_file -- 输出图像文件路径(可选)
    threshold_h -- 杂合度阈值(默认0.5)
    threshold_d -- 深度比例偏差阈值(默认3)
    max_points -- 图中显示的最大点数，超过此数量将进行随机抽样

    返回:
    matplotlib图形对象
    """
    # 对大型数据集进行抽样
    if len(hdplot_table) > max_points:
        logger.info(
            f"数据点数量({len(hdplot_table)})超过最大显示数量({max_points})，进行随机抽样"
        )
        hdplot_table = hdplot_table.sample(max_points, random_state=42)

    # 设置绘图风格
    sns.set(style="whitegrid")

    # 创建图形
    plt.figure(figsize=(12, 10))

    # 绘制散点图
    scatter = plt.scatter(
        hdplot_table["H"],
        hdplot_table["D"],
        alpha=0.6,
        s=15,
        c=hdplot_table["num_hets"],
        cmap="viridis",
    )

    # 添加颜色条
    cbar = plt.colorbar(scatter)
    cbar.set_label("Number of heterozygotes")

    # 添加阈值线
    plt.axhline(y=threshold_d, color="red", linestyle="--", alpha=0.7)
    plt.axhline(y=-threshold_d, color="red", linestyle="--", alpha=0.7)
    plt.axvline(x=threshold_h, color="red", linestyle="--", alpha=0.7)

    # 添加标签和标题
    plt.xlabel("Heterozygosity (H)")
    plt.ylabel("Depth ratio bias (D)")
    plt.title("HDPlot: Identification of paralogous SNPs")

    # 标记可能的旁系同源基因
    potential_paralogs = (hdplot_table["H"] > threshold_h) | (
        abs(hdplot_table["D"]) > threshold_d
    )
    plt.text(
        0.05,
        0.95,
        f"Potential paralogs in sample: {sum(potential_paralogs)}/{len(hdplot_table)}",
        transform=plt.gca().transAxes,
        fontsize=12,
        bbox=dict(facecolor="white", alpha=0.7),
    )

    # 保存图像
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"HDPlot图像已保存至: {output_file}")

    return plt


def filter_paralogs(hdplot_results, threshold_h=0.5, threshold_d=3, output_file=None):
    """
    根据HDPlot结果筛选潜在的旁系同源SNPs

    参数:
    hdplot_results -- HDPlot函数返回的DataFrame
    threshold_h -- 杂合度阈值(默认0.5)
    threshold_d -- 深度比例偏差阈值(默认3)
    output_file -- 输出文件路径(可选)

    返回:
    包含潜在旁系同源SNPs的DataFrame
    """
    # 筛选潜在的旁系同源SNPs
    paralogs = hdplot_results[
        (hdplot_results["H"] > threshold_h) & (abs(hdplot_results["D"]) > threshold_d)
    ]

    logger.info(
        f"共发现 {len(paralogs)} 个潜在的旁系同源SNPs (总数: {len(hdplot_results)})"
    )

    # 保存结果
    if output_file:
        paralogs.to_csv(output_file, index=False)
        logger.info(f"潜在旁系同源SNPs列表已保存至: {output_file}")

    return paralogs


def hdplot_and_filter(vcf_file: Path, output_prefix: Path, chunk_size: int = 10000):
    hdplot_results = hdplot_large_vcf(str(vcf_file), chunk_size, str(output_prefix))
    plot_hdplot(hdplot_results, f"{output_prefix}_hdplot.png")
    paralogs = filter_paralogs(
        hdplot_results, output_file=f"{output_prefix}_paralogs.csv"
    )
    # 输出结果统计
    print(f"分析完成! 共处理 {len(hdplot_results)} 个SNPs")
    print(
        f"发现 {len(paralogs)} 个潜在的旁系同源SNPs ({len(paralogs)/len(hdplot_results)*100:.2f}%)"
    )


# 使用示例
if __name__ == "__main__":
    # 示例参数
    # vcf_file = "large_genome.vcf.gz"  # 大型VCF文件
    # output_prefix = "large_genome"
    # chunk_size = 5000  # 每次处理5000个变异位点
    typer.run(hdplot_and_filter)
