from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import pandas as pd
import pandera as pa
import typer
from attrs import define, field
from pandera.typing import DataFrame, Series


class MutantDBSchema(pa.DataFrameModel):
    chrom: Series[str] = pa.Field(nullable=False)
    pos: Series[int] = pa.Field(nullable=False)
    refer: Series[str] = pa.Field(nullable=False)
    alt: Series[str] = pa.Field(nullable=False)
    type: Series[str] = pa.Field(nullable=False)
    impact: Series[str] = pa.Field(nullable=False)
    gene: Series[str] = pa.Field(nullable=False)
    exon_rank: Series[str] = pa.Field(nullable=True)
    cds_pos: Series[str] = pa.Field(nullable=False)
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


class geneSchema(pa.DataFrameModel):
    gene_id: Series[str] = pa.Field(nullable=False)


class DbDataProcessor:
    def __init__(self):
        self.schema_map = {
            "mutant_db": MutantDBSchema,
            "cds": cdsSchema,
            "gene": geneSchema,
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


def variant_type(ref: str, alt: str) -> str:
    if len(ref) == len(alt):
        return "snp"
    if len(ref) > len(alt):
        return "del"
    return "ins"


def is_ems_snp(ref: str, alt: str) -> bool:
    return f"{ref}{alt}" in ["GA", "CT"]


@define
class DbGeneralStatsProcessor:

    df: pd.DataFrame
    cds_df: pd.DataFrame
    gene_df: pd.DataFrame
    site_df: pd.DataFrame = field(init=False)
    snp_df: pd.DataFrame = field(init=False)

    def __attrs_post_init__(self):
        self.df["variant_type"] = self.df.apply(
            lambda x: variant_type(x["refer"], x["alt"]), axis=1
        )
        self.site_df = self.get_site_df()
        self.snp_df = self.site_df[self.site_df["variant_type"] == "snp"].copy()
        self.snp_df["is_ems_snp"] = self.df.apply(
            lambda x: is_ems_snp(x["refer"], x["alt"]), axis=1
        )

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
        print(gene_mutant_df)
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


def main(db_table: Path, gene_list: Path, cds_bed: Path):
    processor = DbDataProcessor()
    db_df = pd.read_table(db_table)
    gene_df = pd.read_table(gene_list, header=None, names=["gene_id"])
    cds_df = pd.read_table(cds_bed, header=None, names=["chrom", "start", "end"])
    valid_db_df = processor.load_data(db_df, "mutant_db")
    valid_gene_df = processor.load_data(gene_df, "gene")
    valid_cds_df = processor.load_data(cds_df, "cds")
    general_stats_processor = DbGeneralStatsProcessor(
        df=valid_db_df,
        gene_df=valid_gene_df,
        cds_df=valid_cds_df,
    )
    print(general_stats_processor.all_stats())


if __name__ == "__main__":
    typer.run(main)
