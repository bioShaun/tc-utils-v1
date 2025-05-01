from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import pandas as pd
import pandera as pa
import typer
from attrs import define, field
from loguru import logger
from pandera.typing import DataFrame, Series

VA_CHANGE_ORDER = [
    "C:G>T:A",
    "T:A>C:G",
    "Total transitions",
    "C:G>G:C",
    "C:G>A:T",
    "T:A>G:C",
    "T:A>A:T",
    "Total transversions",
    "Transitions/transversions",
]

OUT_VA_TYPES = [
    "stop_gained/lost",
    "start_lost",
    "frameshift_variant",
    "splice_site_variant",
    "missense_variant",
    "synonymous_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
]


class MutantDBSchema(pa.DataFrameModel):
    chrom: Series[str] = pa.Field(nullable=False)
    pos: Series[int] = pa.Field(nullable=False)
    refer: Series[str] = pa.Field(nullable=False)
    alt: Series[str] = pa.Field(nullable=False)
    type: Series[str] = pa.Field(nullable=False)
    impact: Series[str] = pa.Field(nullable=False)
    gene: Series[str] = pa.Field(nullable=False)
    exon_rank: Series[str] = pa.Field(nullable=True)
    cds_pos: Series[str] = pa.Field(nullable=True)
    protein_pos: Series[str] = pa.Field(nullable=True)
    variant: Series[str] = pa.Field(nullable=False)
    sample_id: Series[str] = pa.Field(nullable=False)
    genotype: Series[str] = pa.Field(nullable=False)
    ref_depth: Series[int] = pa.Field(nullable=False)
    alt_depth: Series[int] = pa.Field(nullable=False)

    @classmethod
    @pa.dataframe_check
    def validate_fold_change(cls, df: pd.DataFrame) -> bool:
        """检查 fold change 不为0"""
        return not (df["fold_change"] == 0).any()


class cdsSchema(pa.DataFrameModel):
    chrom: Series[str] = pa.Field(nullable=False)
    start: Series[int] = pa.Field(nullable=False)
    end: Series[int] = pa.Field(nullable=False)


class chrSizeSchema(pa.DataFrameModel):
    chrom: Series[str] = pa.Field(nullable=False)
    chrom_size: Series[int] = pa.Field(nullable=False)


class geneSchema(pa.DataFrameModel):
    gene_id: Series[str] = pa.Field(nullable=False)


class DbDataProcessor:
    def __init__(self):
        self.schema_map = {
            "mutant_db": MutantDBSchema,
            "cds": cdsSchema,
            "gene": geneSchema,
            "chr_size": chrSizeSchema,
        }

    def load_data(self, df: pd.DataFrame, schema_name: str) -> pd.DataFrame:
        """加载并验证数据"""
        schema = self.schema_map.get(schema_name)
        if schema is None:
            raise ValueError(f"无效的schema_name: {schema_name}")
        try:
            # 读取数据
            # df = pd.read_table(file_path)
            # 验证数据
            validated_df = schema.validate(df)
            return validated_df  # type: ignore
        except pa.errors.SchemaError as e:  # type: ignore
            print(f"数据验证失败: {str(e)}")
            raise
        except Exception as e:
            print(f"数据加载失败: {str(e)}")
            raise


# TODO HC LC基因分开统计
@dataclass
class GeneralStats:
    avg_cov: int = 0
    total_variants: int = 0
    total_snp: int = 0
    total_ins: int = 0
    total_del: int = 0
    total_cds_snp: int = 0
    total_cds_ems_snp: int = 0
    total_cds_ems_snp_rate: float = 0
    total_gene_snp: int = 0
    total_gene_ems_snp: int = 0
    total_gene_ems_snp_rate: float = 0
    avg_mutation_per_cds_kb: float = 0
    gene_with_impacted_mutation: int = 0
    total_genes: int = 0
    gene_with_impacted_mutation_rate: float = 0

    def __str__(self) -> str:
        return (
            f"Average coverage: {self.avg_cov}\n"
            f"Total variants: {self.total_variants}\n"
            f"Total SNPs: {self.total_snp}\n"
            f"Total insertions: {self.total_ins}\n"
            f"Total deletions: {self.total_del}\n"
            f"Total CDS SNPs: {self.total_cds_snp}\n"
            f"Total CDS EMS SNPs: {self.total_cds_ems_snp}\n"
            f"Total CDS EMS SNP rate: {self.total_cds_ems_snp_rate}\n"
            f"Total gene SNPs: {self.total_gene_snp}\n"
            f"Total gene EMS SNPs: {self.total_gene_ems_snp}\n"
            f"Total gene EMS SNP rate: {self.total_gene_ems_snp_rate}\n"
            f"Average mutation per CDS kb: {self.avg_mutation_per_cds_kb}\n"
            f"Genes with impacted mutation: {self.gene_with_impacted_mutation}\n"
            f"Total genes: {self.total_genes}\n"
            f"Genes with impacted mutation rate: {self.gene_with_impacted_mutation_rate}\n"
        )

    def to_df(self):
        df_dict = []
        df_dict.append({"stats": "Average mutant coverage", "value": f"{self.avg_cov}"})
        df_dict.append(
            {"stats": "Total number of mutations ", "value": f"{self.total_variants}"}
        )
        df_dict.append({"stats": "Total SNPs", "value": f"{self.total_snp}"})
        df_dict.append({"stats": "Total insertions", "value": f"{self.total_ins}"})
        df_dict.append({"stats": "Total deletions", "value": f"{self.total_del}"})
        df_dict.append({"stats": "Total CDS SNPs", "value": f"{self.total_cds_snp}"})
        df_dict.append(
            {
                "stats": "Total CDS EMS type SNPs (G-A or C-T)",
                "value": f"{self.total_cds_ems_snp}",
            }
        )
        df_dict.append(
            {
                "stats": "CDS EMS type SNPs (G-A or C-T) Percentage",
                "value": f"{self.total_cds_ems_snp_rate * 100:.2f}%",
            }
        )
        df_dict.append(
            {
                "stats": "Total gene [gene body + upstream + downstream] SNPs",
                "value": f"{self.total_gene_snp}",
            }
        )
        df_dict.append(
            {
                "stats": "Total gene EMS type SNPs (G-A or C-T)",
                "value": f"{self.total_gene_ems_snp}",
            }
        )
        df_dict.append(
            {
                "stats": "Total gene EMS type SNPs (G-A or C-T) Percentage",
                "value": f"{self.total_gene_ems_snp_rate * 100:.2f}%",
            }
        )
        df_dict.append(
            {
                "stats": "Average mutation per CDS kb",
                "value": f"{self.avg_mutation_per_cds_kb:.1f}",
            }
        )
        df_dict.append(
            {
                "stats": "Genes with impacted mutation",
                "value": f"{self.gene_with_impacted_mutation}",
            }
        )
        df_dict.append(
            {
                "stats": "Genes with impacted mutation Percentage",
                "value": f"{self.gene_with_impacted_mutation_rate * 100:.2f}%",
            }
        )
        return pd.DataFrame(df_dict)


def variant_type(ref: str, alt: str) -> str:
    if len(ref) == len(alt):
        return "snp"
    if len(ref) > len(alt):
        return "del"
    return "ins"


def is_ems_snp(ref: str, alt: str) -> bool:
    return f"{ref}{alt}" in ["GA", "CT"]


def merge_stop_gain_lost(va_type: str) -> str:
    if va_type in ["stop_gained", "stop_lost"]:
        return "stop_gained/lost"
    if "splice_" in va_type:
        return "splice_site_variant"
    return va_type


def va_change_type(refer_alt):
    """
    将单碱基突变转换为配对突变表示

    参数:
        refer_alt: str, 格式如 'C->T' 或 'G->A'

    返回:
        str: 配对突变表示，如 'C:G>T:A'
    """
    # 定义碱基互补对
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}

    # 解析输入
    ref, alt = refer_alt.split("->")

    # 获取互补碱基
    ref_comp = complement[ref]
    alt_comp = complement[alt]

    # 标准化表示（总是以C>T或T>C或C>G或C>A开头）
    if ref in ["C", "T"]:
        return f"{ref}:{ref_comp}>{alt}:{alt_comp}"
    else:  # ref 是 G 或 A
        return f"{ref_comp}:{ref}>{alt_comp}:{alt}"


@define
class DbGeneralStatsProcessor:

    df: pd.DataFrame
    cds_df: pd.DataFrame
    gene_df: pd.DataFrame
    chr_df: pd.DataFrame
    site_df: pd.DataFrame = field(init=False)
    snp_df: pd.DataFrame = field(init=False)

    def __attrs_post_init__(self):
        self.df["variant_type"] = self.df.apply(
            lambda x: variant_type(x["refer"], x["alt"]), axis=1
        )
        self.df["first_type"] = (
            self.df["type"].map(lambda x: x.split("&")[0]).map(merge_stop_gain_lost)
        )
        self.site_df = self.get_site_df()
        self.snp_df = self.site_df[self.site_df["variant_type"] == "snp"].copy()
        self.snp_df["is_ems_snp"] = self.df.apply(
            lambda x: is_ems_snp(x["refer"], x["alt"]), axis=1
        )
        self.snp_df["refer"] = self.snp_df["refer"].map(lambda x: x[0])
        self.snp_df["alt"] = self.snp_df["alt"].map(lambda x: x[0])
        self.snp_df["refer_to_alt"] = self.snp_df.apply(
            lambda x: f"{x['refer']}->{x['alt']}", axis=1
        )
        self.snp_df["va_change_type"] = self.snp_df["refer_to_alt"].map(va_change_type)

    def get_site_df(self):
        return self.df.drop_duplicates(["chrom", "pos", "refer", "alt"])

    def cal_average_cov(self):
        self.df["total_depth"] = self.df["ref_depth"] + self.df["alt_depth"]
        return int(self.df["total_depth"].mean())

    def cal_total_variants(self):
        return len(self.site_df)

    def cal_total_snp_indel(self) -> Tuple[int, int, int]:
        snp_indel_count = self.site_df["variant_type"].value_counts()
        return (
            snp_indel_count.get("snp", 0),
            snp_indel_count.get("ins", 0),
            snp_indel_count.get("del", 0),
        )

    def cds_variants_stats(self) -> Tuple[int, int, float]:
        cds_snp_df = self.snp_df[~self.snp_df["protein_pos"].isna()]
        ems_snp_count = cds_snp_df["is_ems_snp"].sum()
        return len(cds_snp_df), ems_snp_count, ems_snp_count / len(cds_snp_df)

    def gene_variants_stats(self) -> Tuple[int, int, float]:
        gene_snp_df = self.snp_df[self.snp_df["type"] != "intergenic_region"]
        ems_snp_count = gene_snp_df["is_ems_snp"].sum()
        return len(gene_snp_df), ems_snp_count, ems_snp_count / len(gene_snp_df)

    def impacted_gene_stats(self) -> Tuple[int, float]:
        gene_mutant_df = self.site_df[
            self.site_df["gene"].isin(self.gene_df["gene_id"])
        ]
        impacted_mutant_df = gene_mutant_df[
            gene_mutant_df["impact"].isin(["HIGH", "MODERATE"])
        ]
        impacted_gene_count = len(impacted_mutant_df["gene"].unique())
        return impacted_gene_count, impacted_gene_count / len(self.gene_df)

    def cal_avg_mutation_per_cds_kb(self):
        cds_length = (self.cds_df["end"] - self.cds_df["start"]).sum()
        return len(self.site_df) / (cds_length / 1000)

    def all_stats(self) -> GeneralStats:
        total_variants = self.cal_total_variants()
        total_snp, total_ins, total_del = self.cal_total_snp_indel()
        total_cds_snp, total_cds_ems_snp, total_cds_ems_snp_rate = (
            self.cds_variants_stats()
        )
        total_gene_snp, total_gene_ems_snp, total_gene_ems_snp_rate = (
            self.gene_variants_stats()
        )
        avg_cov = self.cal_average_cov()
        avg_mutation_per_cds_kb = self.cal_avg_mutation_per_cds_kb()
        gene_with_impacted_mutation, gene_with_impacted_mutation_rate = (
            self.impacted_gene_stats()
        )
        return GeneralStats(
            avg_cov=avg_cov,
            total_variants=total_variants,
            total_snp=total_snp,
            total_ins=total_ins,
            total_del=total_del,
            total_cds_snp=total_cds_snp,
            total_cds_ems_snp=total_cds_ems_snp,
            total_cds_ems_snp_rate=total_cds_ems_snp_rate,
            total_gene_snp=total_gene_snp,
            total_gene_ems_snp=total_gene_ems_snp,
            total_gene_ems_snp_rate=total_gene_ems_snp_rate,
            avg_mutation_per_cds_kb=avg_mutation_per_cds_kb,
            gene_with_impacted_mutation=gene_with_impacted_mutation,
            total_genes=len(self.gene_df),
            gene_with_impacted_mutation_rate=gene_with_impacted_mutation_rate,
        )

    def va_count_by_chrom(self) -> pd.DataFrame:
        va_by_chrom = (
            self.site_df.groupby(["chrom", "variant_type"])
            .size()
            .unstack(1)
            .fillna(0)
            .astype(int)
        )
        va_by_chrom["total"] = va_by_chrom.sum(axis=1)
        out_va_by_chrom = va_by_chrom[["snp", "ins", "del", "total"]]
        out_va_by_chrom = self.chr_df.merge(out_va_by_chrom.reset_index())
        out_va_by_chrom.columns = out_va_by_chrom.columns.str.upper()
        return out_va_by_chrom

    def ts_tv_stats(self) -> pd.DataFrame:
        va_change_type_count = self.snp_df["va_change_type"].value_counts()

        transitions = va_change_type_count["C:G>T:A"] + va_change_type_count["T:A>C:G"]
        tansversions = va_change_type_count.sum() - transitions
        ts_tv_ratio = transitions / tansversions

        va_change_type_count_dict = va_change_type_count.to_dict()

        va_change_type_count_dict["Total transitions"] = transitions
        va_change_type_count_dict["Total transversions"] = tansversions
        va_change_type_count_dict["Transitions/transversions"] = f"{ts_tv_ratio:.1f}"

        return (
            pd.DataFrame(va_change_type_count_dict, index=[0])
            .T.loc[VA_CHANGE_ORDER]
            .reset_index()
        )

    def gene_va_type_stats(self) -> pd.DataFrame:
        gene_mutant_df = self.site_df[
            self.site_df["gene"].isin(self.gene_df["gene_id"])
        ]
        return (
            gene_mutant_df.groupby(["first_type"])["gene"]
            .unique()
            .map(len)
            .loc[OUT_VA_TYPES]
            .reset_index()
        )

    def sample_va_stats(self) -> pd.DataFrame:
        sample_va_type_count = (
            self.df.groupby(["sample_id", "first_type"])
            .size()
            .unstack(1)
            .fillna(0)
            .astype("int")
        )[OUT_VA_TYPES].reset_index()
        sample_va_count = (
            self.df.groupby("sample_id").size().rename("total_mutations").reset_index()
        )

        cg_ta_va_count = (
            self.snp_df[self.snp_df["va_change_type"] == "C:G>T:A"]
            .groupby("sample_id")
            .size()
            .rename("C:G>T:A")
            .reset_index()
        )

        return sample_va_count.merge(cg_ta_va_count).merge(sample_va_type_count)


def main(
    db_table: Path, gene_list: Path, cds_bed: Path, chr_size: Path, out_excel: Path
):
    processor = DbDataProcessor()
    db_df = pd.read_table(db_table)
    db_df["chrom"] = db_df["chrom"].astype(str)
    gene_df = pd.read_table(gene_list, header=None, names=["gene_id"])
    cds_df = pd.read_table(cds_bed, header=None, names=["chrom", "start", "end"])
    cds_df["chrom"] = cds_df["chrom"].astype(str)
    chr_df = pd.read_table(chr_size, header=None, names=["chrom", "chrom_size"])
    chr_df["chrom"] = chr_df["chrom"].astype(str)

    logger.info("Start to load data")
    valid_db_df = processor.load_data(db_df, "mutant_db")
    valid_gene_df = processor.load_data(gene_df, "gene")
    valid_cds_df = processor.load_data(cds_df, "cds")
    valid_chr_df = processor.load_data(chr_df, "chr_size")

    general_stats_processor = DbGeneralStatsProcessor(
        df=valid_db_df,
        gene_df=valid_gene_df,
        cds_df=valid_cds_df,
        chr_df=valid_chr_df,
    )
    df_list = []
    logger.info("Start to process table1")
    df_list.append(
        {
            "table": general_stats_processor.all_stats().to_df(),
            "sheet_name": "table1",
            "header": False,
        }
    )
    logger.info("Start to process table2")
    df_list.append(
        {
            "table": general_stats_processor.va_count_by_chrom(),
            "sheet_name": "table2",
            "header": True,
        }
    )
    logger.info("Start to process table3")
    df_list.append(
        {
            "table": general_stats_processor.ts_tv_stats(),
            "sheet_name": "table3",
            "header": False,
        }
    )

    logger.info("Start to process table4")
    df_list.append(
        {
            "table": general_stats_processor.gene_va_type_stats(),
            "sheet_name": "table4",
            "header": False,
        }
    )
    logger.info("Start to process table5")
    df_list.append(
        {
            "table": general_stats_processor.sample_va_stats(),
            "sheet_name": "table5",
            "header": True,
        }
    )
    logger.info("Start to write to excel")
    with pd.ExcelWriter(out_excel) as writer:
        for each_df in df_list:
            each_df["table"].to_excel(
                writer,
                sheet_name=each_df["sheet_name"],
                index=False,
                header=each_df["header"],
            )
    logger.info("Done")


if __name__ == "__main__":
    typer.run(main)
