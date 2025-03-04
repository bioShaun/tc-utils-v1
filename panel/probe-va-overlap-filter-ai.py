from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Optional

import delegator
import pandas as pd
import pandera as pa
import typer
from loguru import logger
from pandera.typing import Series
from pydantic import BaseModel, validator
from tqdm import tqdm

# 配置日志

PROBE_COLUMNS = ["chrom", "probe_start", "probe_end", "id"]


app = typer.Typer()


class ProcessingConfig(BaseModel):
    threads: int
    variant_cutoff: int
    indel_cutoff: int

    @validator("threads")
    def validate_threads(cls, v):
        if v < 1:
            raise ValueError("线程数必须大于0")
        return v

    @validator("variant_cutoff", "indel_cutoff")
    def validate_cutoff(cls, v):
        if v < 0:
            raise ValueError("阈值不能为负数")
        return v


class MutantDBSchema(pa.DataFrameModel):
    chrom: Series[str] = pa.Field(nullable=False)
    probe_start: Series[int] = pa.Field(nullable=False)
    probe_end: Series[int] = pa.Field(nullable=False)
    id: Series[str] = pa.Field(nullable=False)


class TempFileManager:
    def __init__(self):
        self.temp_files = []

    def add(self, file_path: Path):
        self.temp_files.append(file_path)

    def cleanup(self):
        for file_path in self.temp_files:
            try:
                if file_path.exists():
                    file_path.unlink()
            except Exception as e:
                logger.warning(f"清理临时文件失败 {file_path}: {str(e)}")


@dataclass
class AnnDataProcessor:
    ann_table: Path
    validated_df: pd.DataFrame = field(init=False)

    def load_data(self, schema: pa.DataFrameModel) -> None:
        """加载并验证数据"""
        try:
            df = pd.read_table(
                self.ann_table,
            )
            self.validated_df = schema.validate(df)  # type: ignore
            logger.info("数据加载和验证成功")
        except pa.errors.SchemaError as e:  # type: ignore
            logger.error(f"数据验证失败: {str(e)}")
            raise
        except Exception as e:
            logger.error(f"数据加载失败: {str(e)}")
            raise

    def get_validated_data(self) -> pd.DataFrame:
        """缓存验证后的数据"""
        return self.validated_df

    def to_probe_bed(self):
        """将ann_table转换为probe_bed"""
        probe_bed = self.ann_table.with_suffix(".probe.bed")
        self.validated_df.sort_values(["chrom", "probe_start"], inplace=True)
        self.validated_df.to_csv(
            probe_bed, sep="\t", index=False, header=False, columns=PROBE_COLUMNS
        )
        logger.info(f"生成探针BED文件: {probe_bed}")
        return probe_bed


@dataclass
class VcfProcessor:
    vcf_path: Path
    out_dir: Path
    threads: int
    indel_vcf_path: Path = field(init=False)

    def __post_init__(self):
        self.indel_vcf_path = (
            self.out_dir / self.vcf_path.with_suffix(".indel.vcf.gz").name
        )
        self.vcf_bed_path = self.out_dir / self.vcf_path.with_suffix(".bed").name
        self.indel_bed_path = (
            self.out_dir / self.vcf_path.with_suffix(".indel.bed").name
        )

    def indel_filter(self) -> None:
        try:
            if not self.indel_vcf_path.exists():
                cmd = f"bcftools view --exclude-types snps {self.vcf_path} -Oz -o {self.indel_vcf_path} --threads {self.threads}"
                logger.info(cmd)
                result = delegator.run(cmd)
            logger.info("INDEL过滤完成")
        except Exception as e:
            logger.error(f"INDEL过滤失败: {str(e)}")
            raise

    @staticmethod
    def vcf2bed(vcf_path: Path, bed_path: Path) -> None:
        try:
            chunk_size = 10000
            dfs = []

            for chunk in pd.read_table(
                vcf_path,
                header=None,
                names=["chrom", "end"],
                usecols=[0, 1],
                comment="#",
                chunksize=chunk_size,
            ):
                chunk["start"] = chunk["end"] - 1
                chunk["marker"] = 1
                dfs.append(chunk)

            vcf_pos_df = pd.concat(dfs, ignore_index=True)
            vcf_pos_df.to_csv(
                bed_path,
                sep="\t",
                index=False,
                header=False,
                columns=["chrom", "start", "end", "marker"],
            )
            logger.info(f"VCF转换为BED完成: {bed_path}")
        except Exception as e:
            logger.error(f"VCF转换BED失败: {str(e)}")
            raise

    def process_files(self):
        """并行处理多个文件转换任务"""
        self.indel_filter()
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            futures = []
            if not self.indel_bed_path.exists():
                futures.append(
                    executor.submit(
                        self.vcf2bed, self.indel_vcf_path, self.indel_bed_path
                    )
                )
            if not self.vcf_bed_path.exists():
                futures.append(
                    executor.submit(self.vcf2bed, self.vcf_path, self.vcf_bed_path)
                )

            if futures:
                for future in futures:
                    future.result()


@dataclass
class SepVcfProcessor:
    snp_vcf_path: Path
    out_dir: Path
    threads: int
    indel_vcf_path: Optional[Path] = None
    indel_bed_path: Optional[Path] = field(init=False)
    snp_bed_path: Path = field(init=False)

    def __post_init__(self):
        self.snp_bed_path = (
            self.out_dir / self.snp_vcf_path.with_suffix(".snp.bed").name
        )
        if self.indel_vcf_path is not None:
            self.indel_bed_path = (
                self.out_dir / self.indel_vcf_path.with_suffix(".indel.bed").name
            )
            self.all_bed_path = self.snp_bed_path
        else:
            self.indel_bed_path = None
            self.all_bed_path = self.snp_bed_path.with_name("all.bed")


    def merge_bed(self) -> None:
        try:
            if self.indel_bed_path is not None:
                cmd = f"cat {self.snp_bed_path} {self.indel_bed_path} > {self.all_bed_path}"
                logger.info(cmd)
                result = delegator.run(cmd)
            logger.info("BED文件合并完成")
        except Exception as e:
            logger.error(f"合并BED文件失败: {str(e)}")
            raise

    def process_files(self):
        """并行处理多个文件转换任务"""
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            futures = []
            if self.indel_vcf_path is not None:
                if not self.indel_bed_path.exists():
                    futures.append(
                        executor.submit(
                            VcfProcessor.vcf2bed,
                            self.indel_vcf_path,
                            self.indel_bed_path,
                        )
                    )
            if not self.snp_bed_path.exists():
                futures.append(
                    executor.submit(
                        VcfProcessor.vcf2bed, self.snp_vcf_path, self.snp_bed_path
                    )
                )

            if futures:
                for future in futures:
                    future.result()
                self.merge_bed()


def map_variant_to_probe(probe_bed: Path, vcf_bed: Path, col_name: str) -> pd.DataFrame:
    overlap_bed = vcf_bed.with_suffix(".probe.overlap.bed")
    try:
        cmd = f"bedtools map -a {probe_bed} -b {vcf_bed} -c 4 -o sum > {overlap_bed}"
        logger.info(cmd)
        result = delegator.run(cmd)
        if result.return_code != 0:
            raise RuntimeError(f"bedtools 执行失败: {result.err}")

        chunk_size = 10000
        overlap_dfs = []

        total_lines = sum(1 for _ in open(overlap_bed))
        with tqdm(total=total_lines, desc=f"处理{col_name}映射") as pbar:
            for chunk in pd.read_table(
                overlap_bed,
                header=None,
                usecols=[3, 4],
                names=["id", col_name],
                chunksize=chunk_size,
            ):
                chunk[col_name] = chunk[col_name].replace(".", 0).astype(int)
                overlap_dfs.append(chunk)
                pbar.update(len(chunk))

        return pd.concat(overlap_dfs, ignore_index=True)
    except Exception as e:
        logger.error(f"变异映射失败: {str(e)}")
        raise
    # finally:
    #    if overlap_bed.exists():
    #        overlap_bed.unlink()

@app.command()
def seperate_vcf(
    ann_table: Path,
    snp_vcf: Path,
    out_table: Path,
    threads: int = 16,
    variant_cutoff: int = 3,
    indel_cutoff: int = 0,
    indel_vcf: Optional[Path] = None,
    id_list: Optional[Path] = None,
) -> None:
    """主处理流程"""
    out_dir = out_table.parent

    # 验证输入参数
    config = ProcessingConfig(
        threads=threads, variant_cutoff=variant_cutoff, indel_cutoff=indel_cutoff
    )

    if not ann_table.exists():
        raise FileNotFoundError(f"注释文件不存在: {ann_table}")
    if not snp_vcf.exists():
        raise FileNotFoundError(f"VCF文件不存在: {snp_vcf}")

    temp_manager = TempFileManager()

    try:
        # 1. 处理注释数据
        logger.info("开始处理注释数据...")
        ann_data_processor = AnnDataProcessor(ann_table)
        ann_data_processor.load_data(MutantDBSchema)
        probe_bed = ann_data_processor.to_probe_bed()
        temp_manager.add(probe_bed)

        # 2. 处理VCF文件
        logger.info("开始处理VCF文件...")
        vcf_processor = SepVcfProcessor(
            snp_vcf, out_dir, threads=config.threads, indel_vcf_path=indel_vcf
        )
        vcf_processor.process_files()

        # 3. 映射变异到探针
        logger.info("开始变异映射...")
        if vcf_processor.indel_bed_path is not None:
            va_overlap_df = map_variant_to_probe(
                probe_bed, vcf_processor.indel_bed_path, "indel_overlap"
            )

        probe_overlap_df = map_variant_to_probe(
            probe_bed, vcf_processor., "variant_overlap"
        )

        # 4. 合并结果
        if vcf_processor.indel_bed_path is not None:
            add_overlap_df = (
                ann_data_processor.get_validated_data()
                .merge(va_overlap_df)
                .merge(probe_overlap_df)
            )
        else:
            add_overlap_df = ann_data_processor.get_validated_data().merge(probe_overlap_df)
            add_overlap_df["indel_overlap"] = 0
        # 5. 过滤结果
        va_filter = add_overlap_df["variant_overlap"] <= config.variant_cutoff

        indel_filter = add_overlap_df["indel_overlap"] <= config.indel_cutoff
        filter_df = add_overlap_df[va_filter & indel_filter]

        # 6. ID列表过滤(如果提供)
        if id_list and id_list.exists():
            id_set = set(pd.read_table(id_list, header=None)[0])
            filter_df = filter_df[filter_df["id"].isin(id_set)]

        # 7. 保存结果
        filter_df.to_csv(out_table, sep="\t", index=False)
        logger.info(f"处理完成,结果保存至: {out_table}")

    except Exception as e:
        logger.error(f"处理失败: {str(e)}")
        raise
    finally:
        temp_manager.cleanup()

@app.command()
def merged_vcf(
    ann_table: Path,
    vcf: Path,
    out_table: Path,
    threads: int = 16,
    variant_cutoff: int = 3,
    indel_cutoff: int = 0,
    id_list: Optional[Path] = None,
) -> None:
    """主处理流程"""
    out_dir = out_table.parent

    # 验证输入参数
    config = ProcessingConfig(
        threads=threads, variant_cutoff=variant_cutoff, indel_cutoff=indel_cutoff
    )

    if not ann_table.exists():
        raise FileNotFoundError(f"注释文件不存在: {ann_table}")
    if not vcf.exists():
        raise FileNotFoundError(f"VCF文件不存在: {vcf}")

    temp_manager = TempFileManager()

    try:
        # 1. 处理注释数据
        logger.info("开始处理注释数据...")
        ann_data_processor = AnnDataProcessor(ann_table)
        ann_data_processor.load_data(MutantDBSchema)
        probe_bed = ann_data_processor.to_probe_bed()
        temp_manager.add(probe_bed)

        # 2. 处理VCF文件
        logger.info("开始处理VCF文件...")
        vcf_processor = VcfProcessor(vcf, out_dir, threads=config.threads)
        vcf_processor.process_files()

        # 3. 映射变异到探针
        logger.info("开始变异映射...")
        va_overlap_df = map_variant_to_probe(
            probe_bed, vcf_processor.indel_bed_path, "indel_overlap"
        )

        probe_overlap_df = map_variant_to_probe(
            probe_bed, vcf_processor.vcf_bed_path, "variant_overlap"
        )

        # 4. 合并结果
        add_overlap_df = (
            ann_data_processor.get_validated_data()
            .merge(va_overlap_df)
            .merge(probe_overlap_df)
        )

        # 5. 过滤结果
        va_filter = add_overlap_df["variant_overlap"] <= config.variant_cutoff

        indel_filter = add_overlap_df["indel_overlap"] <= config.indel_cutoff
        filter_df = add_overlap_df[va_filter & indel_filter]

        # 6. ID列表过滤(如果提供)
        if id_list and id_list.exists():
            id_set = set(pd.read_table(id_list, header=None)[0])
            filter_df = filter_df[filter_df["id"].isin(id_set)]

        # 7. 保存结果
        filter_df.to_csv(out_table, sep="\t", index=False)
        logger.info(f"处理完成,结果保存至: {out_table}")

    except Exception as e:
        logger.error(f"处理失败: {str(e)}")
        raise
    finally:
        temp_manager.cleanup()


if __name__ == "__main__":
    app()