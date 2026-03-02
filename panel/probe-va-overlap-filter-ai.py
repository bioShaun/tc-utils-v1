from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from inspect import cleandoc
from pathlib import Path
from shlex import quote
from typing import Optional

import delegator
import pandas as pd
import pandera as pa
import typer
from loguru import logger
from pandera.typing import Series
from pydantic import BaseModel, field_validator

try:
    from typing import Annotated
except ImportError:  # pragma: no cover - Python < 3.9 fallback
    from typing_extensions import Annotated

PROBE_COLUMNS = ["chrom", "probe_start", "probe_end", "id"]

MODULE_HELP = cleandoc(
    """
    探针与变异覆盖过滤工具。

    \b
    使用示例:
      python probe-va-overlap-filter-ai.py seperate-vcf ann.tsv snp.vcf.gz out.tsv --threads 8
      python probe-va-overlap-filter-ai.py merged-vcf ann.tsv all.vcf.gz out.tsv --variant-cutoff 3

    \b
    输出格式:
      - out.tsv: 保留注释表原始列，并新增 `variant_overlap` 和 `indel_overlap` 列
    """
)

SEPARATE_VCF_HELP = cleandoc(
    """
    分别输入 SNP VCF 和可选 INDEL VCF，计算探针覆盖并按阈值过滤。

    \b
    使用示例:
      python probe-va-overlap-filter-ai.py seperate-vcf ann.tsv snp.vcf.gz out.tsv --indel-vcf indel.vcf.gz

    \b
    输出格式:
      - out.tsv: 保留输入注释表全部列，并追加 `variant_overlap`、`indel_overlap`
    """
)

MERGED_VCF_HELP = cleandoc(
    """
    输入包含 SNP+INDEL 的混合 VCF，计算探针覆盖并按阈值过滤。

    \b
    使用示例:
      python probe-va-overlap-filter-ai.py merged-vcf ann.tsv all.vcf.gz out.tsv --threads 16

    \b
    输出格式:
      - out.tsv: 保留输入注释表全部列，并追加 `variant_overlap`、`indel_overlap`
    """
)

app = typer.Typer(help=MODULE_HELP, no_args_is_help=True)


def run_command(cmd: str, context: str) -> None:
    """执行外部命令并在失败时抛出包含上下文的异常。"""
    logger.info(cmd)
    result = delegator.run(cmd)
    if result.return_code != 0:
        raise RuntimeError(f"{context}; return_code={result.return_code}; stderr={result.err}")


class ProcessingConfig(BaseModel):
    threads: int
    variant_cutoff: int
    indel_cutoff: int

    @field_validator("threads")
    @classmethod
    def validate_threads(cls, v: int) -> int:
        if v < 1:
            raise ValueError("线程数必须大于0")
        return v

    @field_validator("variant_cutoff", "indel_cutoff")
    @classmethod
    def validate_cutoff(cls, v: int) -> int:
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

    def iter_validated_chunks(
        self, schema: pa.DataFrameModel, chunk_size: int = 100000
    ):
        """按块读取并验证注释表，避免一次性载入超大表格。"""
        try:
            for chunk in pd.read_table(self.ann_table, chunksize=chunk_size):
                chunk["chrom"] = chunk["chrom"].astype(str)
                yield schema.validate(chunk)  # type: ignore
        except pa.errors.SchemaError as e:  # type: ignore
            logger.error(f"数据验证失败: {str(e)}; ann_table={self.ann_table}")
            raise
        except Exception as e:
            logger.error(f"数据加载失败: {str(e)}; ann_table={self.ann_table}")
            raise

    def to_probe_bed(self, schema: pa.DataFrameModel) -> Path:
        """将 ann_table 转换为排序后的 probe BED。"""
        probe_bed = self.ann_table.with_suffix(".probe.bed")
        unsorted_probe_bed = self.ann_table.with_suffix(".probe.unsorted.bed")
        has_rows = False
        try:
            if unsorted_probe_bed.exists():
                unsorted_probe_bed.unlink()

            with unsorted_probe_bed.open("w", encoding="utf-8") as out_handle:
                for chunk in self.iter_validated_chunks(schema=schema):
                    has_rows = True
                    chunk.to_csv(
                        out_handle,
                        sep="\t",
                        index=False,
                        header=False,
                        columns=PROBE_COLUMNS,
                    )

            if not has_rows:
                probe_bed.write_text("", encoding="utf-8")
                return probe_bed

            sort_cmd = (
                f"sort -k1,1 -k2,2n {quote(str(unsorted_probe_bed))} > {quote(str(probe_bed))}"
            )
            run_command(
                sort_cmd,
                f"排序probe BED失败; input={unsorted_probe_bed}; output={probe_bed}",
            )
            logger.info(f"生成探针BED文件: {probe_bed}")
            return probe_bed
        finally:
            if unsorted_probe_bed.exists():
                unsorted_probe_bed.unlink()


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
                cmd = (
                    f"bcftools view --exclude-types snps {quote(str(self.vcf_path))} "
                    f"-Oz -o {quote(str(self.indel_vcf_path))} --threads {self.threads}"
                )
                run_command(
                    cmd,
                    (
                        "INDEL过滤失败"
                        f"; vcf_path={self.vcf_path}; out_path={self.indel_vcf_path}"
                    ),
                )
            logger.info("INDEL过滤完成")
        except Exception as e:
            logger.error(f"INDEL过滤失败: {str(e)}")
            raise

    @staticmethod
    def vcf2bed(vcf_path: Path, bed_path: Path) -> None:
        try:
            chunk_size = 100000
            temp_bed_path = bed_path.with_suffix(".unsorted.bed")
            if temp_bed_path.exists():
                temp_bed_path.unlink()

            with temp_bed_path.open("w", encoding="utf-8") as out_handle:
                for chunk in pd.read_table(
                    vcf_path,
                    header=None,
                    names=["chrom", "end"],
                    usecols=[0, 1],
                    comment="#",
                    chunksize=chunk_size,
                ):
                    chunk["chrom"] = chunk["chrom"].astype(str)
                    chunk["start"] = chunk["end"] - 1
                    chunk["marker"] = 1
                    chunk.to_csv(
                        out_handle,
                        sep="\t",
                        index=False,
                        header=False,
                        columns=["chrom", "start", "end", "marker"],
                    )

            if temp_bed_path.stat().st_size == 0:
                bed_path.write_text("", encoding="utf-8")
            else:
                sort_cmd = (
                    f"sort -k1,1 -k2,2n {quote(str(temp_bed_path))} > {quote(str(bed_path))}"
                )
                run_command(
                    sort_cmd,
                    f"BED排序失败; input={temp_bed_path}; output={bed_path}",
                )
            logger.info(f"VCF转换为BED完成: {bed_path}")
        except Exception as e:
            logger.error(f"VCF转换BED失败: {str(e)}; vcf_path={vcf_path}; bed_path={bed_path}")
            raise
        finally:
            temp_bed_path = bed_path.with_suffix(".unsorted.bed")
            if temp_bed_path.exists():
                temp_bed_path.unlink()

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
    all_bed_path: Path = field(init=False)

    def __post_init__(self):
        self.snp_bed_path = (
            self.out_dir / self.snp_vcf_path.with_suffix(".snp.bed").name
        )
        if self.indel_vcf_path is not None:
            self.indel_bed_path = (
                self.out_dir / self.indel_vcf_path.with_suffix(".indel.bed").name
            )
            self.all_bed_path = self.snp_bed_path.with_name("all.bed")
        else:
            self.indel_bed_path = None
            self.all_bed_path = self.snp_bed_path

    def merge_bed(self) -> None:
        try:
            if self.indel_bed_path is not None:
                cmd = (
                    f"cat {quote(str(self.snp_bed_path))} {quote(str(self.indel_bed_path))} "
                    f"| sort -k1,1 -k2,2n > {quote(str(self.all_bed_path))}"
                )
                run_command(
                    cmd,
                    (
                        "合并BED文件失败"
                        f"; snp_bed={self.snp_bed_path}; indel_bed={self.indel_bed_path}; "
                        f"all_bed={self.all_bed_path}"
                    ),
                )
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
    processed_overlap = overlap_bed.with_suffix(f".{col_name}.tsv")
    try:
        cmd = (
            f"bedtools map -a {quote(str(probe_bed))} -b {quote(str(vcf_bed))} "
            f"-c 4 -o sum > {quote(str(overlap_bed))}"
        )
        run_command(
            cmd,
            f"bedtools 执行失败; probe_bed={probe_bed}; vcf_bed={vcf_bed}",
        )

        has_rows = False
        with processed_overlap.open("w", encoding="utf-8") as out_handle:
            is_first_chunk = True
            for chunk in pd.read_table(
                overlap_bed,
                header=None,
                usecols=[3, 4],
                names=["id", col_name],
                chunksize=100000,
            ):
                has_rows = True
                chunk[col_name] = chunk[col_name].replace(".", 0).astype(int)
                chunk.to_csv(
                    out_handle,
                    sep="\t",
                    index=False,
                    header=is_first_chunk,
                )
                is_first_chunk = False

        if not has_rows:
            return pd.DataFrame(columns=["id", col_name])
        return pd.read_table(processed_overlap, dtype={"id": str, col_name: int})
    except Exception as e:
        logger.error(
            f"变异映射失败: {str(e)}; probe_bed={probe_bed}; vcf_bed={vcf_bed}; col_name={col_name}"
        )
        raise
    finally:
        if overlap_bed.exists():
            overlap_bed.unlink()
        if processed_overlap.exists():
            processed_overlap.unlink()


def write_filtered_output_in_chunks(
    ann_data_processor: AnnDataProcessor,
    out_table: Path,
    config: ProcessingConfig,
    variant_overlap_df: pd.DataFrame,
    indel_overlap_df: Optional[pd.DataFrame] = None,
    id_list: Optional[Path] = None,
) -> None:
    """按块合并与过滤，避免一次性处理超大注释表。"""
    variant_overlap_df = variant_overlap_df.copy()
    variant_overlap_df["id"] = variant_overlap_df["id"].astype(str)
    if indel_overlap_df is not None:
        indel_overlap_df = indel_overlap_df.copy()
        indel_overlap_df["id"] = indel_overlap_df["id"].astype(str)

    base_columns = list(pd.read_table(ann_data_processor.ann_table, nrows=0).columns)
    output_columns = [*base_columns, "indel_overlap", "variant_overlap"]
    write_header = True
    id_set = None
    if id_list and id_list.exists():
        id_set = set(pd.read_table(id_list, header=None)[0].astype(str))

    for ann_chunk in ann_data_processor.iter_validated_chunks(schema=MutantDBSchema):
        ann_chunk["id"] = ann_chunk["id"].astype(str)
        add_overlap_df = ann_chunk.merge(variant_overlap_df, on="id", how="left")

        if indel_overlap_df is not None:
            add_overlap_df = add_overlap_df.merge(indel_overlap_df, on="id", how="left")
        else:
            add_overlap_df["indel_overlap"] = 0

        add_overlap_df["variant_overlap"] = add_overlap_df["variant_overlap"].fillna(0).astype(int)
        add_overlap_df["indel_overlap"] = add_overlap_df["indel_overlap"].fillna(0).astype(int)
        add_overlap_df = add_overlap_df[output_columns]

        va_filter = add_overlap_df["variant_overlap"] <= config.variant_cutoff
        indel_filter = add_overlap_df["indel_overlap"] <= config.indel_cutoff
        filter_df = add_overlap_df[va_filter & indel_filter]

        if id_set is not None:
            filter_df = filter_df[filter_df["id"].astype(str).isin(id_set)]

        if filter_df.empty:
            continue

        filter_df.to_csv(
            out_table,
            sep="\t",
            index=False,
            header=write_header,
            mode="w" if write_header else "a",
        )
        write_header = False

    if write_header:
        pd.DataFrame(columns=output_columns).to_csv(out_table, sep="\t", index=False)


def run_seperate_vcf_pipeline(
    ann_table: Path,
    snp_vcf: Path,
    out_table: Path,
    threads: int = 16,
    variant_cutoff: int = 3,
    indel_cutoff: int = 0,
    indel_vcf: Optional[Path] = None,
    id_list: Optional[Path] = None,
) -> None:
    """执行分离 SNP/INDEL 输入模式的核心处理流程。"""
    out_dir = out_table.parent
    out_dir.mkdir(parents=True, exist_ok=True)

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
        probe_bed = ann_data_processor.to_probe_bed(MutantDBSchema)
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
            probe_bed, vcf_processor.all_bed_path, "variant_overlap"
        )

        # 4. 分块合并与过滤结果，避免一次性处理超大注释表
        write_filtered_output_in_chunks(
            ann_data_processor=ann_data_processor,
            out_table=out_table,
            config=config,
            variant_overlap_df=probe_overlap_df,
            indel_overlap_df=va_overlap_df if vcf_processor.indel_bed_path is not None else None,
            id_list=id_list,
        )
        logger.info(f"处理完成,结果保存至: {out_table}")

    except Exception as e:
        logger.error(f"处理失败: {str(e)}")
        raise
    finally:
        temp_manager.cleanup()


def run_merged_vcf_pipeline(
    ann_table: Path,
    vcf: Path,
    out_table: Path,
    threads: int = 16,
    variant_cutoff: int = 3,
    indel_cutoff: int = 0,
    id_list: Optional[Path] = None,
) -> None:
    """执行混合 VCF 输入模式的核心处理流程。"""
    out_dir = out_table.parent
    out_dir.mkdir(parents=True, exist_ok=True)

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
        probe_bed = ann_data_processor.to_probe_bed(MutantDBSchema)
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

        # 4. 分块合并与过滤结果，避免一次性处理超大注释表
        write_filtered_output_in_chunks(
            ann_data_processor=ann_data_processor,
            out_table=out_table,
            config=config,
            variant_overlap_df=probe_overlap_df,
            indel_overlap_df=va_overlap_df,
            id_list=id_list,
        )
        logger.info(f"处理完成,结果保存至: {out_table}")

    except Exception as e:
        logger.error(f"处理失败: {str(e)}")
        raise
    finally:
        temp_manager.cleanup()


@app.command(help=SEPARATE_VCF_HELP)
def seperate_vcf(
    ann_table: Annotated[Path, typer.Argument(help="注释表文件路径（TSV）")],
    snp_vcf: Annotated[Path, typer.Argument(help="SNP VCF 文件路径（VCF/VCF.GZ）")],
    out_table: Annotated[Path, typer.Argument(help="输出结果文件路径（TSV）")],
    threads: Annotated[int, typer.Option(help="并行线程数，必须大于 0")] = 16,
    variant_cutoff: Annotated[
        int, typer.Option(help="variant_overlap 允许的最大阈值")
    ] = 3,
    indel_cutoff: Annotated[
        int, typer.Option(help="indel_overlap 允许的最大阈值")
    ] = 0,
    indel_vcf: Annotated[
        Optional[Path], typer.Option(help="可选的 INDEL VCF 文件路径")
    ] = None,
    id_list: Annotated[
        Optional[Path], typer.Option(help="可选 ID 白名单文件（单列，无表头）")
    ] = None,
) -> None:
    """CLI 命令入口：分离 SNP/INDEL 输入模式。"""
    try:
        run_seperate_vcf_pipeline(
            ann_table=ann_table,
            snp_vcf=snp_vcf,
            out_table=out_table,
            threads=threads,
            variant_cutoff=variant_cutoff,
            indel_cutoff=indel_cutoff,
            indel_vcf=indel_vcf,
            id_list=id_list,
        )
        typer.echo(f"处理完成: {out_table}")
    except Exception as e:
        logger.exception(f"命令执行失败: {str(e)}")
        typer.echo(f"处理失败: {str(e)}", err=True)
        raise typer.Exit(code=1) from e


@app.command(help=MERGED_VCF_HELP)
def merged_vcf(
    ann_table: Annotated[Path, typer.Argument(help="注释表文件路径（TSV）")],
    vcf: Annotated[Path, typer.Argument(help="混合变异 VCF 文件路径（VCF/VCF.GZ）")],
    out_table: Annotated[Path, typer.Argument(help="输出结果文件路径（TSV）")],
    threads: Annotated[int, typer.Option(help="并行线程数，必须大于 0")] = 16,
    variant_cutoff: Annotated[
        int, typer.Option(help="variant_overlap 允许的最大阈值")
    ] = 3,
    indel_cutoff: Annotated[
        int, typer.Option(help="indel_overlap 允许的最大阈值")
    ] = 0,
    id_list: Annotated[
        Optional[Path], typer.Option(help="可选 ID 白名单文件（单列，无表头）")
    ] = None,
) -> None:
    """CLI 命令入口：混合 VCF 输入模式。"""
    try:
        run_merged_vcf_pipeline(
            ann_table=ann_table,
            vcf=vcf,
            out_table=out_table,
            threads=threads,
            variant_cutoff=variant_cutoff,
            indel_cutoff=indel_cutoff,
            id_list=id_list,
        )
        typer.echo(f"处理完成: {out_table}")
    except Exception as e:
        logger.exception(f"命令执行失败: {str(e)}")
        typer.echo(f"处理失败: {str(e)}", err=True)
        raise typer.Exit(code=1) from e


if __name__ == "__main__":
    app()
