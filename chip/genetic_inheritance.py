import gzip
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd
import typer
from cyvcf2 import VCF
from tqdm import tqdm

VALID_MODES = {"f1", "f2"}
MAX_INCONSISTENT_DETAILS = 100


def normalize_mode(mode: str) -> str:
    """规范化并校验遗传模式参数。"""
    normalized_mode = mode.strip().lower()
    if normalized_mode not in VALID_MODES:
        raise ValueError(f"无效模式: {mode}，仅支持: {', '.join(sorted(VALID_MODES))}")
    return normalized_mode


class VCFGeneticAnalyzer:
    def __init__(self, mode: str = "f1"):
        self.vcf_data = None
        self.family_combinations = None
        self.results = None
        self.sample_names = None
        self.variant_info = None
        self.mode = normalize_mode(mode)

    def parse_vcf_header(self, file_path: str) -> List[str]:
        """解析VCF文件头部，获取样品名称"""
        sample_names = []

        try:
            # 判断是否为压缩文件
            open_func = gzip.open if file_path.endswith(".gz") else open
            mode = "rt" if file_path.endswith(".gz") else "r"

            with open_func(file_path, mode) as f:
                for line in f:
                    if line.startswith("#CHROM"):
                        # VCF头部格式: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample1 sample2 ...
                        headers = line.strip().split("\t")
                        if len(headers) > 9:  # 标准VCF至少有10列
                            sample_names = headers[9:]  # 从第10列开始是样品名
                        break
                    elif not line.startswith("#"):
                        break

            return sample_names

        except Exception as e:
            print(f"❌ 解析VCF头部失败: {str(e)}")
            return []

    def load_vcf_data(
        self, file_path: Path, max_alleles: int = 2, max_variants: Optional[int] = None
    ) -> pd.DataFrame:
        """使用cyvcf2高效读取VCF中的基因型矩阵和变异信息"""
        try:
            print(f"📁 正在加载VCF文件: {file_path} (cyvcf2)")
            vcf = VCF(file_path)
            sample_names = list(vcf.samples)
            self.sample_names = sample_names
            print(f"✓ 发现 {len(sample_names)} 个样品")
            print(
                f"样品名称: {sample_names[:10]}{'...' if len(sample_names) > 10 else ''}"
            )
            records = []
            gt_data = []
            variant_ids = []
            for idx, variant in enumerate(vcf):
                if max_variants is not None and idx >= max_variants:
                    break
                # 生成variant_id，兼容多ALT
                if len(variant.ALT) > max_alleles - 1:
                    print(
                        f"⚠️ 变异位点 {variant.CHROM}:{variant.POS} 有超过 {max_alleles - 1} 个等位基因，跳过"
                    )
                    continue
                alt_str = ",".join(variant.ALT)
                var_id = f"{variant.CHROM}:{variant.POS}:{variant.REF}:{alt_str}"
                variant_ids.append(var_id)
                records.append(
                    [
                        variant.CHROM,
                        variant.POS,
                        variant.ID or ".",
                        variant.REF,
                        alt_str,
                        variant.QUAL,
                        variant.FILTER or ".",
                        variant.INFO,
                        variant.FORMAT,
                    ]
                )
                # 获取标准GT字段， cyvcf2.genotypes: shape=(samples, 3) e.g. [[0, 1, True], ...]
                gt_row = []
                for n, gt in enumerate(variant.genotypes):
                    # gt: e.g. [0, 1, True], or [None, None, False] (缺失)
                    if gt[0] is None or gt[1] is None or gt[0] == -1 or gt[1] == -1:
                        gt_str = "./."
                    # elif gt[0] != gt[1]:
                    #     af = variant.gt_alt_freqs[n]
                    #     if af < 0.2 or af > 0.8:
                    #         gt_str = "./."
                    #     gt_str = f"{gt[0]}/{gt[1]}"
                    else:
                        sep = "/"  # 实际上cyvcf2能区分phase，这里只用'/'即可
                        gt_str = f"{gt[0]}{sep}{gt[1]}"
                    gt_row.append(gt_str)
                gt_data.append(gt_row)
                if idx > 0 and idx % 10000 == 0:
                    print(f"  已加载 {idx} 个变异位点...")
            if not records:
                print("❌ 未找到有效的变异数据")
                return pd.DataFrame()
            # 转为DataFrame
            vcf_columns = [
                "CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
            ]
            variant_info_df = pd.DataFrame(
                records, columns=vcf_columns, index=variant_ids
            )
            self.variant_info = variant_info_df
            gt_df = pd.DataFrame(gt_data, columns=sample_names, index=variant_ids)
            self.vcf_data = gt_df
            print(f"✓ 成功加载 {len(gt_df)} 个变异位点")
            return gt_df
        except Exception as e:
            print(f"❌ 加载VCF文件失败: {str(e)}")
            return pd.DataFrame()

    def load_family_combinations(self, file_path: Path) -> pd.DataFrame:
        """加载家系组合文件"""
        try:
            if file_path.suffix in (".xlsx", ".xls"):
                df = pd.read_excel(file_path)
            elif file_path.suffix == ".tsv":
                df = pd.read_csv(file_path, sep="\t")
            elif file_path.suffix == ".csv":
                df = pd.read_csv(file_path)
            else:
                print("❌ 不支持的文件格式")
                return pd.DataFrame()

            print(f"✓ 成功加载家系组合数据: {len(df)} 个家系")
            print(f"列名: {list(df.columns)}")

            # 标准化列名
            column_mapping = {}
            for col in df.columns:
                col_lower = col.lower()
                if "father" in col_lower or "父" in col_lower or "dad" in col_lower:
                    column_mapping[col] = "father"
                elif "mother" in col_lower or "母" in col_lower or "mom" in col_lower:
                    column_mapping[col] = "mother"
                elif (
                    "child" in col_lower
                    or "子" in col_lower
                    or "offspring" in col_lower
                ):
                    column_mapping[col] = "child"

            if len(column_mapping) < 3:
                if len(df.columns) >= 3:
                    column_mapping = {
                        df.columns[0]: "father",
                        df.columns[1]: "mother",
                        df.columns[2]: "child",
                    }
                    print("⚠️  自动将前三列识别为: father, mother, child")
                else:
                    raise ValueError("无法识别父本-母本-子代列")

            df.rename(columns=column_mapping, inplace=True)

            # 检查样品是否存在于VCF数据中
            if self.sample_names is not None:
                missing_samples = []
                for _, row in df.iterrows():
                    for role in ["father", "mother", "child"]:
                        if row[role] not in self.sample_names:
                            missing_samples.append(f"{row[role]} ({role})")

                if missing_samples:
                    print("⚠️  以下样品在VCF文件中未找到:")
                    for sample in missing_samples[:10]:
                        print(f"    {sample}")
                    if len(missing_samples) > 10:
                        print(f"    ... 还有 {len(missing_samples)-10} 个样品")

            self.family_combinations = df[["father", "mother", "child"]].copy()
            return self.family_combinations

        except Exception as e:
            print(f"❌ 加载家系组合数据失败: {str(e)}")
            return pd.DataFrame()

    def parse_vcf_genotype(self, gt_string: str) -> List[str]:
        """
        解析VCF格式的基因型
        支持格式: 0/0, 0/1, 1/1, 0|1, ./., 等
        """
        if pd.isna(gt_string) or gt_string == "" or gt_string == ".":
            return list()

        # 提取基因型部分（FORMAT字段的第一个值通常是GT）
        gt_parts = str(gt_string).split(":")
        if not gt_parts:
            return list()

        gt_field = gt_parts[0]  # 基因型字段

        # 处理缺失值
        if gt_field == "./." or gt_field == ".|." or gt_field == ".":
            return list()

        # 分离等位基因
        if "/" in gt_field:
            alleles = gt_field.split("/")
        elif "|" in gt_field:
            alleles = gt_field.split("|")
        else:
            # 单个等位基因
            alleles = [gt_field]

        # 过滤掉缺失的等位基因
        valid_alleles = [a for a in alleles if a != "."]

        return valid_alleles

    def resolve_alleles(self, variant_id: str, allele_indices: "Iterable[str]") -> list:
        """
        将等位基因索引转换为实际的核苷酸序列（保留顺序，不去重）
        0 = REF, 1 = ALT1, 2 = ALT2, etc.
        """
        resolved_alleles = []
        try:
            variant_info = self.variant_info.loc[variant_id]
            ref_allele = variant_info["REF"]
            alt_alleles = variant_info["ALT"].split(",")
            for idx in allele_indices:
                if idx == "0":
                    resolved_alleles.append(ref_allele)
                elif idx.isdigit():
                    alt_idx = int(idx) - 1
                    if 0 <= alt_idx < len(alt_alleles):
                        resolved_alleles.append(alt_alleles[alt_idx])
                # else 忽略‘.’或其它
            return resolved_alleles
        except Exception:
            # 如果解析失败，返回原始索引（以列表方式）
            return list(allele_indices)

    def check_variant_consistency(
        self, variant_id: str, father_gt: str, mother_gt: str, child_gt: str
    ) -> Tuple[bool, str]:
        """检查单个位点遗传一致性，支持F1和F2两种模式"""
        try:
            # 解析基因型为等位基因索引
            father_indices = self.parse_vcf_genotype(father_gt)
            mother_indices = self.parse_vcf_genotype(mother_gt)
            child_indices = self.parse_vcf_genotype(child_gt)

            # 检查缺失
            if not father_indices or not mother_indices or not child_indices:
                return False, "缺失基因型数据"

            # 转为实际等位基因（列表，保留顺序）
            father_alleles = list(self.resolve_alleles(variant_id, father_indices))
            mother_alleles = list(self.resolve_alleles(variant_id, mother_indices))
            child_alleles = list(self.resolve_alleles(variant_id, child_indices))

            # 要求每个样本必须有两个等位基因（仅支持二倍体）
            if (
                len(father_alleles) != 2
                or len(mother_alleles) != 2
                or len(child_alleles) != 2
            ):
                return False, "不是二倍体基因型"

            if self.mode == "f1":
                # F1模式: 子代必须从父、母各继承一个等位基因
                possible_children = set()
                for f in father_alleles:
                    for m in mother_alleles:
                        # 不考虑等位基因顺序（纯合、杂合都允许），用frozenset或排序tuple统一表达
                        possible_children.add(tuple(sorted([f, m])))

                # 实际子代基因型
                child_tuple = tuple(sorted(child_alleles))

                if child_tuple in possible_children:
                    return True, "遗传一致"
                else:
                    return (
                        False,
                        f"子代({child_alleles})不可能由父({father_alleles})母({mother_alleles})组合得到",
                    )
            else:
                # F2模式: 子代两个等位基因只需存在于基因池(父∪母)中
                gene_pool = set(father_alleles) | set(mother_alleles)

                child_in_pool = [allele in gene_pool for allele in child_alleles]

                if all(child_in_pool):
                    return True, "遗传一致"
                else:
                    missing_alleles = [a for a, in_pool in zip(child_alleles, child_in_pool) if not in_pool]
                    return (
                        False,
                        f"子代等位基因{missing_alleles}不在基因池({list(gene_pool)})中"
                    )
        except Exception as e:
            return False, f"分析错误: {str(e)}"

    def analyze_family_consistency(self, chunk_size: int = 1000) -> pd.DataFrame:
        """
        分析所有家系的遗传一致性
        chunk_size: 分批处理的变异位点数量
        """
        if self.vcf_data is None or self.family_combinations is None:
            print("❌ 请先加载VCF数据和家系组合数据")
            return pd.DataFrame()

        results = []
        variants = self.vcf_data.index
        total_variants = len(variants)

        print(
            f"🔍 开始分析 {len(self.family_combinations)} 个家系，{total_variants} 个变异位点..."
        )

        for idx, family in self.family_combinations.iterrows():
            father_id = family["father"]
            mother_id = family["mother"]
            child_id = family["child"]

            print(f"  分析家系 {idx+1}: {father_id} × {mother_id} → {child_id}")

            # 检查样品是否存在
            missing_samples = []
            if father_id not in self.vcf_data.columns:
                missing_samples.append(f"父本({father_id})")
            if mother_id not in self.vcf_data.columns:
                missing_samples.append(f"母本({mother_id})")
            if child_id not in self.vcf_data.columns:
                missing_samples.append(f"子代({child_id})")

            if missing_samples:
                results.append(
                    {
                        "family_id": idx + 1,
                        "father": father_id,
                        "mother": mother_id,
                        "child": child_id,
                        "total_variants": total_variants,
                        "consistent_variants": 0,
                        "inconsistent_variants": 0,
                        "missing_data_variants": total_variants,
                        "consistency_rate": 0.0,
                        "status": f"样品缺失: {', '.join(missing_samples)}",
                    }
                )
                continue

            # 分批处理变异位点
            consistent_count = 0
            inconsistent_count = 0
            missing_count = 0
            inconsistent_details = []

            for i in tqdm(range(0, total_variants, chunk_size)):
                end_idx = min(i + chunk_size, total_variants)
                chunk_variants = variants[i:end_idx]

                for variant_id in chunk_variants:
                    father_gt = self.vcf_data.loc[variant_id, father_id]
                    mother_gt = self.vcf_data.loc[variant_id, mother_id]
                    child_gt = self.vcf_data.loc[variant_id, child_id]

                    is_consistent, explanation = self.check_variant_consistency(
                        variant_id, father_gt, mother_gt, child_gt
                    )

                    if "缺失" in explanation:
                        missing_count += 1
                    elif is_consistent:
                        consistent_count += 1
                    else:
                        inconsistent_count += 1
                        if len(inconsistent_details) < MAX_INCONSISTENT_DETAILS:
                            inconsistent_details.append(
                                {
                                    "variant_id": variant_id,
                                    "chrom": self.variant_info.loc[variant_id, "CHROM"],
                                    "pos": self.variant_info.loc[variant_id, "POS"],
                                    "ref": self.variant_info.loc[variant_id, "REF"],
                                    "alt": self.variant_info.loc[variant_id, "ALT"],
                                    "father_gt": father_gt,
                                    "mother_gt": mother_gt,
                                    "child_gt": child_gt,
                                    "issue": explanation,
                                }
                            )
            # 计算一致性比例
            valid_variants = total_variants - missing_count
            consistency_rate = (
                consistent_count / valid_variants if valid_variants > 0 else 0
            )

            family_result = {
                "family_id": idx + 1,
                "father": father_id,
                "mother": mother_id,
                "child": child_id,
                "total_variants": total_variants,
                "consistent_variants": consistent_count,
                "inconsistent_variants": inconsistent_count,
                "missing_data_variants": missing_count,
                "consistency_rate": consistency_rate,
                "status": (
                    "正常"
                    if inconsistent_count == 0
                    else f"{inconsistent_count}个变异位点不一致"
                ),
            }

            if inconsistent_details:
                family_result["inconsistent_details"] = inconsistent_details

            results.append(family_result)

        self.results = pd.DataFrame(results)
        return self.results

    def generate_detailed_report(self, output_dir: str = "."):
        """生成详细报告"""
        if self.results is None:
            print("❌ 请先运行分析")
            return

        # 生成汇总报告
        self.generate_summary_report()

        # 生成详细的变异位点报告
        problem_families = self.results[self.results["inconsistent_variants"] > 0]

        if len(problem_families) > 0:
            print("\n📝 生成详细变异位点报告...")

            for _, family in problem_families.iterrows():
                if "inconsistent_details" in family and family["inconsistent_details"]:
                    family_id = family["family_id"]
                    filename = (
                        f"{output_dir}/family_{family_id}_inconsistent_variants.csv"
                    )

                    details_df = pd.DataFrame(family["inconsistent_details"])
                    details_df.to_csv(filename, index=False)
                    print(f"  家系 {family_id} 不一致变异位点详情: {filename}")

    def generate_summary_report(self):
        """生成汇总报告"""
        if self.results is None:
            print("❌ 请先运行分析")
            return

        print("\n" + "=" * 80)
        print("VCF家系遗传一致性分析报告")
        print("=" * 80)

        total_families = len(self.results)
        fully_consistent = len(self.results[self.results["inconsistent_variants"] == 0])
        has_issues = total_families - fully_consistent

        print(f"总家系数量: {total_families}")
        print(
            f"完全一致家系: {fully_consistent} ({fully_consistent/total_families:.1%})"
        )
        print(f"存在问题家系: {has_issues} ({has_issues/total_families:.1%})")

        if has_issues > 0:
            avg_consistency = self.results["consistency_rate"].mean()
            print(f"平均一致性比例: {avg_consistency:.1%}")

            # 不一致变异位点统计
            total_inconsistent = self.results["inconsistent_variants"].sum()
            total_analyzed = (
                self.results["total_variants"].sum()
                - self.results["missing_data_variants"].sum()
            )
            overall_error_rate = (
                total_inconsistent / total_analyzed if total_analyzed > 0 else 0
            )

            print(f"总体错误率: {overall_error_rate:.3%}")
            print(f"总不一致变异位点: {total_inconsistent}")

            print("\n问题家系Top10:")
            print("-" * 60)

            worst_families = self.results.nsmallest(10, "consistency_rate")
            for _, family in worst_families.iterrows():
                if family["inconsistent_variants"] > 0:
                    print(
                        f"家系 {family['family_id']}: {family['father']} × {family['mother']} → {family['child']}"
                    )
                    print(
                        f"  一致性: {family['consistency_rate']:.1%} "
                        + f"({family['consistent_variants']}/{family['total_variants']})"
                    )
                    print(f"  不一致变异位点: {family['inconsistent_variants']}")
                    print()

    def save_results(self, output_file: Path):
        """保存分析结果"""
        if self.results is None:
            print("❌ 请先运行分析")
            return

        try:
            output_data = self.results.copy()
            if "inconsistent_details" in output_data.columns:
                output_data = output_data.drop("inconsistent_details", axis=1)

            if output_file.suffix in (".xlsx", ".xls"):
                output_data.to_excel(output_file, index=False)
            elif output_file.suffix == ".tsv":
                output_data.to_csv(output_file, index=False, sep="\t")
            elif output_file.suffix == ".csv":
                output_data.to_csv(output_file, index=False)
            else:
                print("❌ 不支持的输出文件格式")
                return

            print(f"✓ 结果已保存到: {output_file}")

        except Exception as e:
            print(f"❌ 保存失败: {str(e)}")


def main(
    vcf_file: Path,
    family_file: Path,
    output_file: Path,
    max_variants: Optional[int] = None,
    max_alleles: int = 2,
    mode: str = typer.Option(
        "f1",
        "--mode",
        help="遗传一致性检查模式: f1=直接亲子(F1代), f2=原始亲本到F2代",
    ),
):
    """主程序"""
    try:
        normalized_mode = normalize_mode(mode)
    except ValueError as exc:
        raise typer.BadParameter(str(exc), param_hint="--mode") from exc

    analyzer = VCFGeneticAnalyzer(mode=normalized_mode)

    print("VCF文件家系遗传一致性分析工具")
    print(f"模式: {'F1(直接亲子)' if normalized_mode == 'f1' else 'F2(原始亲本到F2代)'}")
    print("=" * 60)

    # 输入文件路径
    # 询问是否限制变异位点数量（用于大文件测试）
    if max_variants:
        try:
            max_variants = int(max_variants)
        except ValueError:
            max_variants = None
    else:
        max_variants = None

    # 加载数据
    print("\n📁 加载数据文件...")
    vcf_data = analyzer.load_vcf_data(
        vcf_file, max_alleles=max_alleles, max_variants=max_variants
    )
    family_data = analyzer.load_family_combinations(family_file)

    if vcf_data.empty or family_data.empty:
        print("❌ 数据加载失败，程序退出")
        return

    # 运行分析
    print("\n🧬 开始遗传一致性分析...")
    results = analyzer.analyze_family_consistency()

    if not results.empty:
        # 显示报告
        analyzer.generate_detailed_report()

        # 保存结果
        analyzer.save_results(output_file)

        print("\n✅ 分析完成！")
    else:
        print("❌ 分析失败")


if __name__ == "__main__":
    typer.run(main)
