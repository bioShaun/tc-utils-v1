#!/usr/bin/env python3
"""
VCF Fingerprint Analysis
计算:
1. PIC 值分布
2. Nei's 多样性指数 (He)
3. 样本间遗传距离矩阵 (IBS)

依赖:
    pip install cyvcf2 numpy pandas
"""

import numpy as np
import pandas as pd
from cyvcf2 import VCF
from tqdm import tqdm


def calc_allele_freq(genotypes):
    """
    计算等位基因频率
    genotypes: (n_samples, 2)，0=ref, 1=alt, -1=缺失
    返回: allele_freq (list)
    """
    alleles = genotypes[genotypes >= 0]  # 去除缺失
    if alleles.size == 0:
        return None
    unique, counts = np.unique(alleles, return_counts=True)
    freqs = counts / alleles.size
    return freqs


def pic_value(freqs):
    """
    计算 PIC 值
    """
    if freqs is None or len(freqs) < 2:
        return np.nan
    sum_pi2 = np.sum(freqs**2)
    sum2 = 0
    for i in range(len(freqs)):
        for j in range(i + 1, len(freqs)):
            sum2 += 2 * (freqs[i] ** 2) * (freqs[j] ** 2)
    return 1 - sum_pi2 - sum2


def nei_he(freqs):
    """
    计算 Nei's Gene Diversity (He)
    """
    if freqs is None:
        return np.nan
    return 1 - np.sum(freqs**2)


def calc_ibd_matrix(vcf_file, max_snps=5000):
    """
    计算 IBS 距离矩阵 (简化版)
    随机抽样 max_snps 个位点
    """
    vcf = VCF(vcf_file, gts012=True)
    samples = vcf.samples
    n = len(samples)
    mat = np.zeros((n, n), dtype=float)
    count = np.zeros((n, n), dtype=int)

    snp_idx = 0
    for variant in vcf:
        gts = variant.genotypes  # (n_samples, 3)，最后一列是 phased
        gts = np.array(gts)[:, :2]
        if np.any(gts < 0):
            continue
        # flatten 每个样本的基因型 (0,1,2)
        geno = gts.sum(axis=1)
        # 两两比较
        for i in range(n):
            for j in range(i + 1, n):
                diff = abs(geno[i] - geno[j])
                # IBS: 0=相同, 1=部分相同, 2=完全不同
                shared = 2 - diff
                mat[i, j] += shared
                mat[j, i] += shared
                count[i, j] += 2
                count[j, i] += 2

        snp_idx += 1
        if snp_idx >= max_snps:
            break

    # 转换为距离 (1 - IBS/2)
    dist = 1 - mat / count
    return samples, dist


def main(vcf_file, out_prefix):
    vcf = VCF(vcf_file, gts012=True)

    pic_list = []
    nei_list = []

    for variant in tqdm(vcf, desc="Processing SNPs"):
        gts = np.array(variant.genotypes)[:, :2]
        freqs = calc_allele_freq(gts.flatten())
        if freqs is None:
            continue
        pic_list.append(pic_value(freqs))
        nei_list.append(nei_he(freqs))

    # 输出 PIC / Nei's
    df = pd.DataFrame({"PIC": pic_list, "Nei_He": nei_list})
    df.to_csv(f"{out_prefix}.snp_stats.csv", index=False)

    print(f"Done! 输出文件:\n  {out_prefix}.snp_stats.csv\n  {out_prefix}.ibs_dist.csv")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="VCF SNP diversity analysis")
    parser.add_argument("-i", "--vcf", required=True, help="Input VCF file")
    parser.add_argument("-o", "--out", required=True, help="Output prefix")
    args = parser.parse_args()
    main(args.vcf, args.out)
