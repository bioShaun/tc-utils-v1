"""
SNP位点最小化选择工具
目标: 用最少的SNP位点区分所有个体

支持多种算法:
1. 贪心算法 (快速，适用于大数据集)
2. 暴力搜索 (精确，适用于小数据集)
3. 遗传算法 (平衡，适用于中等数据集)
"""

import os
import random
import sys
import time
import warnings
from collections import defaultdict
from itertools import combinations, product
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

# VCF处理库
try:
    from cyvcf2 import VCF

    CYVCF2_AVAILABLE = True
except ImportError:
    CYVCF2_AVAILABLE = False
    warnings.warn("cyvcf2 not available, install with: pip install cyvcf2")

try:
    import pysam

    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False
    warnings.warn("pysam not available, install with: pip install pysam")

# 如果都没有安装，提供安装提示
if not CYVCF2_AVAILABLE and not PYSAM_AVAILABLE:
    print("请安装VCF处理库:")
    print("pip install cyvcf2 pysam")
    print("或者:")
    print("conda install -c bioconda cyvcf2 pysam")


class SNPMinimizer:
    def __init__(
        self,
        snp_data: np.ndarray,
        individual_ids: List[str] = None,
        snp_ids: List[str] = None,
    ):
        """
        初始化SNP最小化器

        Parameters:
        -----------
        snp_data : np.ndarray
            SNP数据矩阵，行为个体，列为SNP位点
            值通常为0,1,2表示基因型，或其他编码方式
        individual_ids : List[str], optional
            个体ID列表
        snp_ids : List[str], optional
            SNP位点ID列表
        """
        self.snp_data = snp_data
        self.n_individuals, self.n_snps = snp_data.shape

        self.individual_ids = individual_ids or [
            f"Individual_{i}" for i in range(self.n_individuals)
        ]
        self.snp_ids = snp_ids or [f"SNP_{i}" for i in range(self.n_snps)]

        print(f"数据概览: {self.n_individuals}个个体, {self.n_snps}个SNP位点")

    def _can_distinguish_all(self, snp_indices: List[int]) -> bool:
        """检查给定的SNP位点组合是否能区分所有个体"""
        if not snp_indices:
            return False

        # 提取选定SNP位点的基因型
        selected_data = self.snp_data[:, snp_indices]

        # 将每个个体的基因型转换为元组用于比较
        genotypes = [tuple(row) for row in selected_data]

        # 检查是否所有基因型都是唯一的
        return len(set(genotypes)) == self.n_individuals

    def greedy_selection(self) -> Tuple[List[int], List[str]]:
        """
        贪心算法选择最少SNP位点
        每次选择能最大化区分个体数量的SNP位点
        """
        print("运行贪心算法...")
        start_time = time.time()

        selected_snps = []
        remaining_snps = list(range(self.n_snps))

        while not self._can_distinguish_all(selected_snps) and remaining_snps:
            best_snp = None
            max_distinction = 0

            for snp_idx in remaining_snps:
                test_snps = selected_snps + [snp_idx]
                # 计算当前组合能区分的个体对数量
                distinction_count = self._count_distinguished_pairs(test_snps)

                if distinction_count > max_distinction:
                    max_distinction = distinction_count
                    best_snp = snp_idx

            if best_snp is not None:
                selected_snps.append(best_snp)
                remaining_snps.remove(best_snp)
                print(
                    f"选择SNP {best_snp} ({self.snp_ids[best_snp]}), 当前已选择{len(selected_snps)}个"
                )
            else:
                break

        end_time = time.time()
        print(f"贪心算法完成，耗时: {end_time - start_time:.2f}秒")

        selected_names = [self.snp_ids[i] for i in selected_snps]
        return selected_snps, selected_names

    def _count_distinguished_pairs(self, snp_indices: List[int]) -> int:
        """计算给定SNP组合能区分的个体对数量"""
        if not snp_indices:
            return 0

        selected_data = self.snp_data[:, snp_indices]
        distinguished_pairs = 0

        for i in range(self.n_individuals):
            for j in range(i + 1, self.n_individuals):
                if not np.array_equal(selected_data[i], selected_data[j]):
                    distinguished_pairs += 1

        return distinguished_pairs

    def brute_force_search(self, max_snps: int = None) -> Tuple[List[int], List[str]]:
        """
        暴力搜索最小SNP位点集合
        从1个SNP开始，逐步增加直到能区分所有个体
        """
        if max_snps is None:
            max_snps = min(20, self.n_snps)  # 限制搜索范围避免计算爆炸

        print(f"运行暴力搜索算法 (最大搜索{max_snps}个SNP)...")
        start_time = time.time()

        for k in range(1, max_snps + 1):
            print(f"搜索{k}个SNP的组合...")

            total_combinations = self._combination_count(self.n_snps, k)
            if total_combinations > 1000000:  # 如果组合数太大，跳过
                print(f"  {k}个SNP的组合数({total_combinations})太大，跳过")
                continue

            for snp_combination in combinations(range(self.n_snps), k):
                if self._can_distinguish_all(list(snp_combination)):
                    end_time = time.time()
                    print(f"暴力搜索完成，耗时: {end_time - start_time:.2f}秒")

                    selected_names = [self.snp_ids[i] for i in snp_combination]
                    return list(snp_combination), selected_names

        end_time = time.time()
        print(f"暴力搜索未找到解，耗时: {end_time - start_time:.2f}秒")
        return [], []

    def _combination_count(self, n: int, k: int) -> int:
        """计算组合数C(n,k)"""
        if k > n or k < 0:
            return 0
        if k == 0 or k == n:
            return 1

        result = 1
        for i in range(min(k, n - k)):
            result = result * (n - i) // (i + 1)
        return result

    def genetic_algorithm(
        self, population_size: int = 100, generations: int = 100
    ) -> Tuple[List[int], List[str]]:
        """
        遗传算法寻找最优SNP位点组合
        """
        print(f"运行遗传算法 (种群大小: {population_size}, 世代数: {generations})...")
        start_time = time.time()

        # 初始化种群
        population = []
        for _ in range(population_size):
            # 随机选择SNP数量和位点
            num_snps = random.randint(1, min(20, self.n_snps))
            snps = random.sample(range(self.n_snps), num_snps)
            population.append(snps)

        best_solution = None
        best_fitness = float("inf")

        for generation in range(generations):
            # 计算适应度
            fitness_scores = []
            for individual in population:
                fitness = self._fitness_function(individual)
                fitness_scores.append(fitness)

                if (
                    self._can_distinguish_all(individual)
                    and len(individual) < best_fitness
                ):
                    best_fitness = len(individual)
                    best_solution = individual.copy()

            if generation % 20 == 0:
                print(
                    f"  第{generation}代，当前最佳解: {best_fitness if best_solution else '未找到'}"
                )

            # 选择、交叉、变异
            population = self._evolve_population(population, fitness_scores)

        end_time = time.time()
        print(f"遗传算法完成，耗时: {end_time - start_time:.2f}秒")

        if best_solution:
            selected_names = [self.snp_ids[i] for i in best_solution]
            return best_solution, selected_names
        else:
            return [], []

    def _fitness_function(self, snp_indices: List[int]) -> float:
        """
        适应度函数
        优先考虑能区分所有个体的解，其次考虑SNP数量少的解
        """
        if not snp_indices:
            return float("inf")

        if self._can_distinguish_all(snp_indices):
            return len(snp_indices)  # 能完全区分时，SNP越少越好
        else:
            # 不能完全区分时，惩罚分数 = SNP数量 + 未区分个体对数量
            total_pairs = self.n_individuals * (self.n_individuals - 1) // 2
            distinguished_pairs = self._count_distinguished_pairs(snp_indices)
            undistinguished_pairs = total_pairs - distinguished_pairs
            return len(snp_indices) + undistinguished_pairs * 10

    def _evolve_population(
        self, population: List[List[int]], fitness_scores: List[float]
    ) -> List[List[int]]:
        """进化种群：选择、交叉、变异"""
        new_population = []

        # 保留最优个体（精英选择）
        sorted_indices = sorted(
            range(len(fitness_scores)), key=lambda i: fitness_scores[i]
        )
        elite_size = len(population) // 10
        for i in range(elite_size):
            new_population.append(population[sorted_indices[i]].copy())

        # 生成新个体
        while len(new_population) < len(population):
            # 锦标赛选择
            parent1 = self._tournament_selection(population, fitness_scores)
            parent2 = self._tournament_selection(population, fitness_scores)

            # 交叉
            child = self._crossover(parent1, parent2)

            # 变异
            child = self._mutate(child)

            new_population.append(child)

        return new_population

    def _tournament_selection(
        self,
        population: List[List[int]],
        fitness_scores: List[float],
        tournament_size: int = 3,
    ) -> List[int]:
        """锦标赛选择"""
        tournament_indices = random.sample(
            range(len(population)), min(tournament_size, len(population))
        )
        best_idx = min(tournament_indices, key=lambda i: fitness_scores[i])
        return population[best_idx].copy()

    def _crossover(self, parent1: List[int], parent2: List[int]) -> List[int]:
        """交叉操作"""
        all_snps = list(set(parent1 + parent2))
        child_size = random.randint(1, min(len(all_snps), 15))
        return random.sample(all_snps, child_size)

    def _mutate(self, individual: List[int], mutation_rate: float = 0.1) -> List[int]:
        """变异操作"""
        if random.random() < mutation_rate:
            if random.random() < 0.5 and len(individual) > 1:
                # 删除一个SNP
                individual.remove(random.choice(individual))
            else:
                # 添加一个新SNP
                available_snps = [i for i in range(self.n_snps) if i not in individual]
                if available_snps:
                    individual.append(random.choice(available_snps))
        return individual

    def analyze_results(
        self, selected_snps: List[int], selected_names: List[str]
    ) -> None:
        """分析和展示结果"""
        if not selected_snps:
            print("未找到有效解决方案")
            return

        print(f"\n=== 结果分析 ===")
        print(f"选择的SNP数量: {len(selected_snps)}")
        print(f"选择的SNP位点: {selected_names}")
        print(f"选择的SNP索引: {selected_snps}")

        # 验证结果
        if self._can_distinguish_all(selected_snps):
            print("✓ 验证通过：所选SNP位点能够区分所有个体")
        else:
            print("✗ 验证失败：所选SNP位点无法区分所有个体")

        # 展示每个个体的基因型
        print(f"\n各个体在所选SNP位点的基因型:")
        selected_data = self.snp_data[:, selected_snps]

        for i, individual_id in enumerate(self.individual_ids):
            genotype = selected_data[i]
            print(f"{individual_id}: {genotype}")


class VCFParser:
    """高效VCF文件解析器 - 使用cyvcf2和pysam"""

    @staticmethod
    def parse_vcf_cyvcf2(
        vcf_file: str,
        max_snps: int = None,
        min_maf: float = 0.05,
        max_missing: float = 0.2,
    ) -> Tuple[np.ndarray, List[str], List[str]]:
        """
        使用cyvcf2解析VCF文件 (推荐，速度最快)

        Parameters:
        -----------
        vcf_file : str
            VCF文件路径
        max_snps : int, optional
            最大SNP数量限制
        min_maf : float, default=0.05
            最小等位基因频率
        max_missing : float, default=0.2
            最大缺失率 (0.2 = 20%)
        """
        if not CYVCF2_AVAILABLE:
            raise ImportError("cyvcf2 not installed. Install with: pip install cyvcf2")

        print(f"使用cyvcf2解析VCF文件: {vcf_file}")

        if not os.path.exists(vcf_file):
            raise FileNotFoundError(f"VCF文件不存在: {vcf_file}")

        # 打开VCF文件
        vcf = VCF(vcf_file)

        # 获取样本名称
        individual_ids = list(vcf.samples)
        n_individuals = len(individual_ids)
        print(
            f"发现 {n_individuals} 个个体: {individual_ids[:5]}{'...' if n_individuals > 5 else ''}"
        )

        snp_data_list = []
        snp_ids = []
        snp_info = []
        processed_count = 0
        skipped_count = 0

        for variant in vcf:
            if max_snps and processed_count >= max_snps:
                break

            # 跳过非SNP变异
            if not variant.is_snp:
                continue

            # 跳过多等位基因位点
            if len(variant.ALT) > 1:
                skipped_count += 1
                continue

            # 获取基因型数据
            genotypes = variant.gt_types  # 0: hom_ref, 1: het, 2: hom_alt, 3: unknown

            # 计算缺失率
            missing_rate = np.sum(genotypes == 3) / len(genotypes)
            if missing_rate > max_missing:
                skipped_count += 1
                continue

            # 处理缺失数据：替换为最常见的基因型
            if missing_rate > 0:
                valid_genotypes = genotypes[genotypes != 3]
                if len(valid_genotypes) > 0:
                    most_common = np.bincount(valid_genotypes).argmax()
                    genotypes[genotypes == 3] = most_common

            # 计算MAF
            allele_count = np.sum(genotypes)
            total_alleles = len(genotypes) * 2
            maf = (
                min(allele_count / total_alleles, 1 - allele_count / total_alleles)
                if total_alleles > 0
                else 0
            )

            # 过滤稀有变异
            if maf < min_maf:
                skipped_count += 1
                continue

            # 构建SNP ID
            snp_id = variant.ID if variant.ID else f"{variant.CHROM}:{variant.POS}"

            # 存储数据
            snp_data_list.append(genotypes)
            snp_ids.append(snp_id)
            snp_info.append(
                {
                    "chrom": variant.CHROM,
                    "pos": variant.POS,
                    "ref": variant.REF,
                    "alt": variant.ALT[0],
                    "maf": maf,
                    "missing_rate": missing_rate,
                }
            )

            processed_count += 1

            if processed_count % 5000 == 0:
                print(f"已处理 {processed_count} 个SNP位点 (跳过 {skipped_count} 个)")

        vcf.close()

        if not snp_data_list:
            raise ValueError("未找到符合条件的SNP数据")

        # 转换为numpy数组 (行为个体，列为SNP)
        snp_data = np.array(snp_data_list).T

        print(f"cyvcf2解析完成:")
        print(f"  - 有效SNP: {snp_data.shape[1]} 个")
        print(f"  - 跳过SNP: {skipped_count} 个")
        print(f"  - 个体数量: {snp_data.shape[0]} 个")
        print(
            f"  - 基因型分布: {dict(zip(['0/0', '0/1', '1/1'], np.bincount(snp_data.flatten(), minlength=3)))}"
        )
        print(
            f"  - MAF范围: {min([info['maf'] for info in snp_info]):.3f} - {max([info['maf'] for info in snp_info]):.3f}"
        )

        return snp_data, individual_ids, snp_ids

    @staticmethod
    def parse_vcf_pysam(
        vcf_file: str,
        max_snps: int = None,
        min_maf: float = 0.05,
        max_missing: float = 0.2,
    ) -> Tuple[np.ndarray, List[str], List[str]]:
        """
        使用pysam解析VCF文件 (备选方案)

        Parameters:
        -----------
        vcf_file : str
            VCF文件路径
        max_snps : int, optional
            最大SNP数量限制
        min_maf : float, default=0.05
            最小等位基因频率
        max_missing : float, default=0.2
            最大缺失率
        """
        if not PYSAM_AVAILABLE:
            raise ImportError("pysam not installed. Install with: pip install pysam")

        print(f"使用pysam解析VCF文件: {vcf_file}")

        if not os.path.exists(vcf_file):
            raise FileNotFoundError(f"VCF文件不存在: {vcf_file}")

        # 打开VCF文件
        vcf_file_obj = pysam.VariantFile(vcf_file)

        # 获取样本名称
        individual_ids = list(vcf_file_obj.header.samples)
        n_individuals = len(individual_ids)
        print(
            f"发现 {n_individuals} 个个体: {individual_ids[:5]}{'...' if n_individuals > 5 else ''}"
        )

        snp_data_list = []
        snp_ids = []
        snp_info = []
        processed_count = 0
        skipped_count = 0

        for record in vcf_file_obj.fetch():
            if max_snps and processed_count >= max_snps:
                break

            # 跳过非SNP或多等位基因
            if (
                len(record.ref) != 1
                or len(record.alts) != 1
                or len(record.alts[0]) != 1
            ):
                skipped_count += 1
                continue

            # 提取基因型
            genotypes = []
            missing_count = 0

            for sample in individual_ids:
                gt = record.samples[sample]["GT"]
                if None in gt:  # 缺失数据
                    genotypes.append(-1)  # 标记为缺失
                    missing_count += 1
                else:
                    # 转换为加性编码
                    genotypes.append(sum(gt))

            # 检查缺失率
            missing_rate = missing_count / len(genotypes)
            if missing_rate > max_missing:
                skipped_count += 1
                continue

            # 处理缺失数据
            genotypes = np.array(genotypes)
            if missing_count > 0:
                valid_genotypes = genotypes[genotypes != -1]
                if len(valid_genotypes) > 0:
                    most_common = np.bincount(valid_genotypes).argmax()
                    genotypes[genotypes == -1] = most_common

            # 计算MAF
            allele_count = np.sum(genotypes)
            total_alleles = len(genotypes) * 2
            maf = (
                min(allele_count / total_alleles, 1 - allele_count / total_alleles)
                if total_alleles > 0
                else 0
            )

            # 过滤稀有变异
            if maf < min_maf:
                skipped_count += 1
                continue

            # 构建SNP ID
            snp_id = record.id if record.id else f"{record.chrom}:{record.pos}"

            # 存储数据
            snp_data_list.append(genotypes)
            snp_ids.append(snp_id)
            snp_info.append(
                {
                    "chrom": record.chrom,
                    "pos": record.pos,
                    "ref": record.ref,
                    "alt": record.alts[0],
                    "maf": maf,
                    "missing_rate": missing_rate,
                }
            )

            processed_count += 1

            if processed_count % 5000 == 0:
                print(f"已处理 {processed_count} 个SNP位点 (跳过 {skipped_count} 个)")

        vcf_file_obj.close()

        if not snp_data_list:
            raise ValueError("未找到符合条件的SNP数据")

        # 转换为numpy数组 (行为个体，列为SNP)
        snp_data = np.array(snp_data_list).T

        print(f"pysam解析完成:")
        print(f"  - 有效SNP: {snp_data.shape[1]} 个")
        print(f"  - 跳过SNP: {skipped_count} 个")
        print(f"  - 个体数量: {snp_data.shape[0]} 个")
        print(
            f"  - 基因型分布: {dict(zip(['0/0', '0/1', '1/1'], np.bincount(snp_data.flatten(), minlength=3)))}"
        )
        print(
            f"  - MAF范围: {min([info['maf'] for info in snp_info]):.3f} - {max([info['maf'] for info in snp_info]):.3f}"
        )

        return snp_data, individual_ids, snp_ids

    @staticmethod
    def parse_vcf(
        vcf_file: str,
        max_snps: int = None,
        min_maf: float = 0.05,
        max_missing: float = 0.2,
        parser: str = "auto",
    ) -> Tuple[np.ndarray, List[str], List[str]]:
        """
        统一VCF解析接口

        Parameters:
        -----------
        vcf_file : str
            VCF文件路径
        max_snps : int, optional
            最大SNP数量限制
        min_maf : float, default=0.05
            最小等位基因频率
        max_missing : float, default=0.2
            最大缺失率
        parser : str, default='auto'
            解析器选择: 'cyvcf2', 'pysam', 'auto'
        """
        if parser == "auto":
            if CYVCF2_AVAILABLE:
                parser = "cyvcf2"
            elif PYSAM_AVAILABLE:
                parser = "pysam"
            else:
                raise ImportError(
                    "Neither cyvcf2 nor pysam is available. Please install at least one."
                )

        if parser == "cyvcf2":
            return VCFParser.parse_vcf_cyvcf2(vcf_file, max_snps, min_maf, max_missing)
        elif parser == "pysam":
            return VCFParser.parse_vcf_pysam(vcf_file, max_snps, min_maf, max_missing)
        else:
            raise ValueError(
                f"Unknown parser: {parser}. Choose 'cyvcf2', 'pysam', or 'auto'"
            )


def load_vcf_data(
    vcf_file: str,
    max_snps: int = None,
    min_maf: float = 0.05,
    max_missing: float = 0.2,
    parser: str = "auto",
) -> Tuple[np.ndarray, List[str], List[str]]:
    """
    从VCF文件加载SNP数据 (使用专业库)

    Parameters:
    -----------
    vcf_file : str
        VCF文件路径
    max_snps : int, optional
        最大SNP数量限制
    min_maf : float, default=0.05
        最小等位基因频率
    max_missing : float, default=0.2
        最大缺失率
    parser : str, default='auto'
        解析器: 'cyvcf2' (推荐), 'pysam', 'auto'
    """
    return VCFParser.parse_vcf(vcf_file, max_snps, min_maf, max_missing, parser)


def load_example_data() -> Tuple[np.ndarray, List[str], List[str]]:
    """生成示例数据用于测试"""
    np.random.seed(42)

    n_individuals = 10
    n_snps = 50

    # 生成随机SNP数据 (0, 1, 2 表示基因型)
    snp_data = np.random.choice(
        [0, 1, 2], size=(n_individuals, n_snps), p=[0.4, 0.4, 0.2]
    )

    # 确保每个个体都有独特的基因型（至少在某些位点上）
    for i in range(1, n_individuals):
        # 随机修改几个位点确保个体间有差异
        diff_positions = np.random.choice(
            n_snps, size=np.random.randint(1, 5), replace=False
        )
        for pos in diff_positions:
            snp_data[i, pos] = (snp_data[i, pos] + np.random.randint(1, 3)) % 3

    individual_ids = [f"Sample_{i+1}" for i in range(n_individuals)]
    snp_ids = [f"rs{1000000+i}" for i in range(n_snps)]

    return snp_data, individual_ids, snp_ids


def create_example_vcf(filename: str = "example.vcf") -> str:
    """创建示例VCF文件用于测试"""
    vcf_content = """##fileformat=VCFv4.2
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=242193529>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2	Sample3	Sample4	Sample5	Sample6	Sample7	Sample8	Sample9	Sample10
1	1000	rs1000	A	G	60	PASS	AF=0.3	GT	0/0	0/1	1/1	0/0	0/1	0/0	1/1	0/1	0/0	0/1
1	2000	rs2000	C	T	60	PASS	AF=0.4	GT	0/1	1/1	0/0	0/1	1/1	0/1	0/0	1/1	0/1	0/0
1	3000	rs3000	G	A	60	PASS	AF=0.5	GT	1/1	0/0	0/1	1/1	0/0	1/1	0/1	0/0	1/1	0/1
2	1000	rs4000	T	C	60	PASS	AF=0.2	GT	0/0	0/0	0/1	0/0	0/0	1/1	0/0	0/1	0/0	0/0
2	2000	rs5000	A	T	60	PASS	AF=0.6	GT	0/1	1/1	1/1	0/1	1/1	0/0	1/1	1/1	0/1	1/1
2	3000	rs6000	C	G	60	PASS	AF=0.35	GT	0/0	0/1	0/0	1/1	0/1	0/0	0/1	1/1	0/1	0/0
1	4000	rs7000	G	T	60	PASS	AF=0.25	GT	0/0	0/0	1/1	0/1	0/0	0/0	0/1	0/0	1/1	0/1
1	5000	rs8000	A	C	60	PASS	AF=0.45	GT	1/1	0/1	0/0	1/1	0/1	1/1	0/0	0/1	0/0	1/1
2	4000	rs9000	T	G	60	PASS	AF=0.15	GT	0/0	0/1	0/0	0/0	0/0	0/1	0/0	0/0	1/1	0/0
2	5000	rs10000	C	A	60	PASS	AF=0.55	GT	0/1	1/1	0/1	1/1	1/1	0/0	1/1	0/1	1/1	0/1
"""

    with open(filename, "w") as f:
        f.write(vcf_content)

    print(f"示例VCF文件已创建: {filename}")
    return filename


def main():
    """主函数示例"""
    print("SNP位点最小化选择工具 - 支持VCF输入")
    print("=" * 50)

    # 检查命令行参数
    if len(sys.argv) > 1:
        vcf_file = sys.argv[1]
        max_snps = int(sys.argv[2]) if len(sys.argv) > 2 else None
        min_maf = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05

        if not os.path.exists(vcf_file):
            print(f"错误：VCF文件不存在: {vcf_file}")
            return

        print(f"从VCF文件加载数据: {vcf_file}")
        try:
            snp_data, individual_ids, snp_ids = load_vcf_data(
                vcf_file, max_snps, min_maf
            )
        except Exception as e:
            print(f"VCF文件解析错误: {e}")
            return
    else:
        # 创建并使用示例VCF文件
        print("未提供VCF文件，使用示例数据")
        example_vcf = create_example_vcf()
        snp_data, individual_ids, snp_ids = load_vcf_data(
            example_vcf, max_snps=50, min_maf=0.1
        )

    # 创建SNP最小化器
    minimizer = SNPMinimizer(snp_data, individual_ids, snp_ids)

    # 选择合适的算法
    if minimizer.n_snps <= 20:
        # 小数据集：使用暴力搜索获得最优解
        print("\n数据集较小，使用暴力搜索算法获得最优解:")
        selected_snps, selected_names = minimizer.brute_force_search()
        minimizer.analyze_results(selected_snps, selected_names)

    elif minimizer.n_snps <= 1000:
        # 中等数据集：使用遗传算法
        print("\n使用遗传算法:")
        selected_snps, selected_names = minimizer.genetic_algorithm(
            population_size=100, generations=100
        )
        minimizer.analyze_results(selected_snps, selected_names)

        # 也运行贪心算法进行比较
        print("\n" + "=" * 50)
        print("贪心算法结果 (用于比较):")
        greedy_snps, greedy_names = minimizer.greedy_selection()
        minimizer.analyze_results(greedy_snps, greedy_names)

    else:
        # 大数据集：仅使用贪心算法
        print("\n数据集较大，使用贪心算法:")
        selected_snps, selected_names = minimizer.greedy_selection()
        minimizer.analyze_results(selected_snps, selected_names)

    # 清理示例文件
    if len(sys.argv) == 1 and os.path.exists("example.vcf"):
        os.remove("example.vcf")
        print("\n示例文件已清理")


def main_with_vcf(
    vcf_file: str,
    max_snps: int = None,
    min_maf: float = 0.05,
    max_missing: float = 0.2,
    algorithm: str = "auto",
    parser: str = "auto",
):
    """
    直接使用VCF文件的主函数 (使用专业库)

    Parameters:
    -----------
    vcf_file : str
        VCF文件路径
    max_snps : int, optional
        最大SNP数量限制
    min_maf : float, default=0.05
        最小等位基因频率
    max_missing : float, default=0.2
        最大缺失率
    algorithm : str, default='auto'
        算法选择: 'greedy', 'brute_force', 'genetic', 'auto'
    parser : str, default='auto'
        VCF解析器: 'cyvcf2', 'pysam', 'auto'
    """
    print(f"分析VCF文件: {vcf_file}")
    print(f"使用解析器: {parser}")

    # 加载VCF数据
    try:
        snp_data, individual_ids, snp_ids = load_vcf_data(
            vcf_file, max_snps, min_maf, max_missing, parser
        )
    except Exception as e:
        print(f"VCF文件加载失败: {e}")
        return None, None

    # 创建分析器
    minimizer = SNPMinimizer(snp_data, individual_ids, snp_ids)

    # 选择算法
    if algorithm == "auto":
        if minimizer.n_snps <= 20:
            algorithm = "brute_force"
        elif minimizer.n_snps <= 1000:
            algorithm = "genetic"
        else:
            algorithm = "greedy"

    print(f"使用算法: {algorithm}")

    # 运行算法
    if algorithm == "greedy":
        selected_snps, selected_names = minimizer.greedy_selection()
    elif algorithm == "brute_force":
        selected_snps, selected_names = minimizer.brute_force_search()
    elif algorithm == "genetic":
        selected_snps, selected_names = minimizer.genetic_algorithm()
    else:
        raise ValueError(f"未知算法: {algorithm}")

    # 分析结果
    minimizer.analyze_results(selected_snps, selected_names)

    return selected_snps, selected_names


if __name__ == "__main__":
    main()

# 使用说明:
"""
=== VCF文件输入支持 ===

1. 命令行使用:
   python script.py your_file.vcf [max_snps] [min_maf]
   
   示例:
   python script.py data.vcf 1000 0.05
   python script.py data.vcf.gz  # 支持压缩文件

2. 编程接口使用:
   # 方法1: 直接调用
   selected_snps, selected_names = main_with_vcf('your_file.vcf')
   
   # 方法2: 分步骤
   snp_data, individual_ids, snp_ids = load_vcf_data('your_file.vcf')
   minimizer = SNPMinimizer(snp_data, individual_ids, snp_ids)
   selected_snps, selected_names = minimizer.greedy_selection()

3. VCF文件要求:
   - 标准VCF v4.x格式
   - 支持压缩文件 (.vcf.gz)
   - 需要GT (genotype) 字段
   - 自动过滤多等位基因位点
   - 自动过滤稀有变异 (MAF < 0.05)
   - 自动处理缺失数据

4. 基因型编码:
   - 0/0 → 0 (纯合子参考)
   - 0/1, 1/0 → 1 (杂合子)
   - 1/1 → 2 (纯合子变异)

5. 算法选择建议:
   - SNP < 20: 暴力搜索 (最优解)
   - SNP < 1000: 遗传算法 (平衡)
   - SNP >= 1000: 贪心算法 (快速)

6. 输出结果:
   - 最小SNP位点集合
   - 每个个体的基因型
   - 结果验证信息

=== 专业工具推荐 ===
对于大规模数据，也可考虑:
1. PLINK: plink --vcf input.vcf --indep-pairwise 50 10 0.8
2. TASSEL: 内建SNP筛选功能
3. VCFtools: vcftools --vcf input.vcf --maf 0.05 --max-maf 0.95
"""
