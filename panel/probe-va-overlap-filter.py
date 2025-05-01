from dataclasses import dataclass, field
from pathlib import Path

import delegator
import pandas as pd
import pandera as pa
import typer
from pandera.typing import Series

PROBE_COLUMNS = ["chrom", "probe_start", "probe_end", "id"]


class MutantDBSchema(pa.DataFrameModel):
    chrom: Series[str] = pa.Field(nullable=False)
    probe_start: Series[int] = pa.Field(nullable=False)
    probe_end: Series[int] = pa.Field(nullable=False)
    id: Series[str] = pa.Field(nullable=False)


@dataclass
class AnnDataProcessor:

    ann_table: Path
    validated_df: pd.DataFrame = field(init=False)

    def load_data(self, schema: pa.DataFrameModel) -> None:
        """加载并验证数据"""
        try:
            # 读取数据
            df = pd.read_table(self.ann_table)
            # 验证数据
            self.validated_df = schema.validate(df)
        except pa.errors.SchemaError as e:  # type: ignore
            print(f"数据验证失败: {str(e)}")
            raise
        except Exception as e:
            print(f"数据加载失败: {str(e)}")
            raise

    def to_probe_bed(self):
        """将ann_table转换为probe_bed"""
        probe_bed = self.ann_table.with_suffix(".probe.bed")
        self.validated_df.sort_values(["chrom", "probe_start"], inplace=True)
        self.validated_df.to_csv(
            probe_bed, sep="\t", index=False, header=False, columns=PROBE_COLUMNS
        )
        return probe_bed


@dataclass
class VcfProcessor:
    vcf_path: Path
    threads: int
    indel_vcf_path: Path = field(init=False)

    def __post_init__(self):
        self.indel_vcf_path = self.vcf_path.with_suffix(".indel.vcf.gz")
        self.vcf_bed_path = self.vcf_path.with_suffix(".bed")
        self.indel_bed_path = self.vcf_path.with_suffix(".indel.bed")

    def indel_filter(self) -> None:
        delegator.run(
            f"bcftools --exclude-types snps {self.vcf_path} -Oz -o {self.indel_vcf_path} --threads {self.threads}"
        )

    @staticmethod
    def vcf2bed(vcf_path: Path, bed_path: Path) -> None:
        vcf_pos_df = pd.read_table(
            vcf_path, header=None, names=["chrom", "end"], usecols=[0, 1], comment="#"
        )
        vcf_pos_df["start"] = vcf_pos_df["end"] - 1
        vcf_pos_df["marker"] = 1
        vcf_pos_df.to_csv(
            bed_path,
            sep="\t",
            index=False,
            header=False,
            columns=["chrom", "start", "end", "marker"],
        )


def map_variant_to_probe(probe_bed: Path, vcf_bed: Path, col_name: str) -> pd.DataFrame:
    overlap_bed = vcf_bed.with_suffix(".probe.overlap.bed")
    delegator.run(
        f"bedtools map -a {probe_bed} -b {vcf_bed} -c 4 -o sum > {overlap_bed}"
    )
    overlap_df = pd.read_table(
        overlap_bed, header=None, usecols=[3, 4], names=["id", col_name]
    )
    overlap_df[col_name] = overlap_df[col_name].replace(".", 0).astype(int)
    return overlap_df


def main(
    ann_table: Path,
    vcf: Path,
    out_table: Path,
    threads: int = 16,
    variant_cutoff: int = 3,
    indel_cutoff: int = 0,
    id_list: Path = typer.Option(None),
) -> None:

    ann_data_processor = AnnDataProcessor(ann_table)
    ann_data_processor.load_data(MutantDBSchema)
    probe_bed = ann_data_processor.to_probe_bed()

    vcf_processor = VcfProcessor(vcf, threads=threads)
    vcf_processor.indel_filter()
    vcf_processor.vcf2bed(vcf_processor.indel_vcf_path, vcf_processor.indel_bed_path)
    vcf_processor.vcf2bed(vcf_processor.vcf_path, vcf_processor.vcf_bed_path)

    va_overlap_df = map_variant_to_probe(
        probe_bed, vcf_processor.indel_bed_path, "variant_overlap"
    )

    probe_overlap_df = map_variant_to_probe(
        probe_bed, vcf_processor.vcf_bed_path, "indel_overlap"
    )

    add_overlap_df = ann_data_processor.validated_df.merge(va_overlap_df).merge(
        probe_overlap_df
    )

    va_filter = add_overlap_df["variant_overlap"] <= variant_cutoff
    indel_filter = add_overlap_df["indel_overlap"] <= indel_cutoff

    filter_df = add_overlap_df[va_filter & indel_filter]
    if id_list:
        filter_df = filter_df[filter_df["id"].isin(id_list.read_text().split("\n"))]
    filter_df.to_csv(out_table, sep="\t", index=False)
