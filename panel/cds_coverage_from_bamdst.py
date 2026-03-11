"""
根据 bamdst 的 region.tsv.gz 输出，评估 CDS 区域的覆盖度。

使用示例:
    # 默认输出
    python panel/cds_coverage_from_bamdst.py cds_cov_dir/ output.tsv

    # 带 split bed 的版本
    python panel/cds_coverage_from_bamdst.py cds_cov_dir/ output.tsv --split-bed split.bed

    # 输出 xlsx 格式
    python panel/cds_coverage_from_bamdst.py cds_cov_dir/ output.tsv --xlsx

输出格式:
    - TSV 或 XLSX 文件，包含：chrom, pos(end), min_cov, max_cov, mean_cov, quantile_cov, coverage_{cov}x
"""

from inspect import cleandoc
from pathlib import Path
from typing import Annotated

import pandas as pd
import typer
from loguru import logger

MODULE_HELP = cleandoc(
    """
    根据 bamdst 的 region.tsv.gz 输出，评估 CDS 区域的覆盖度。

    \b
    使用示例:
      1) 默认输出 TSV
         python panel/cds_coverage_from_bamdst.py cds_cov_dir/ output.tsv

    \b
      2) 带 split bed 的版本
         python panel/cds_coverage_from_bamdst.py cds_cov_dir/ output.tsv --split-bed split.bed

    \b
      3) 同时输出 xlsx 格式
         python panel/cds_coverage_from_bamdst.py cds_cov_dir/ output.tsv --xlsx

    \b
    输出格式:
      - chrom: 染色体
      - pos(end): 区域结束位置（原 end 列）
      - min_cov: 最小覆盖度
      - max_cov: 最大覆盖度
      - mean_cov: 平均覆盖度
      - quantile_cov: 中位数覆盖度
      - coverage_{{cov}}x: 各阈值下的覆盖比例
    """
)

app = typer.Typer(help=MODULE_HELP, no_args_is_help=True)


def merge_chr(df: pd.DataFrame, split_bed: Path) -> pd.DataFrame:
    if not split_bed.exists():
        logger.error(f"split bed 文件不存在: {split_bed}")
        raise typer.Exit(code=1)

    split_bed_df = pd.read_csv(
        split_bed,
        header=None,
        names=["new_chrom", "offset", "offset_end", "chrom"],
        sep="\t",
    )
    merged_df = df.merge(split_bed_df, on="chrom", how="left", validate="many_to_one")
    unmatched = merged_df["offset"].isna().sum()
    if unmatched > 0:
        logger.warning(f"有 {unmatched} 个区域未匹配到 split bed 中的染色体")
    merged_df["new_start"] = merged_df["start"] + merged_df["offset"]
    merged_df["new_end"] = merged_df["end"] + merged_df["offset"]
    merged_df = merged_df.drop(["chrom", "offset", "offset_end", "start", "end"], axis=1)
    merged_df = merged_df.rename(
        columns={"new_chrom": "chrom", "new_start": "start", "new_end": "end"},
    )
    return merged_df


def load_bed_files(bed_dir: Path, cov_cutoff: float | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not bed_dir.exists():
        logger.error(f"输入目录不存在: {bed_dir}")
        raise typer.Exit(code=1)

    bed_list = sorted(bed_dir.glob("*/region.tsv.gz"))
    if not bed_list:
        logger.error(f"在 {bed_dir} 中未找到 */region.tsv.gz 文件")
        raise typer.Exit(code=1)

    df_list: list[pd.DataFrame] = []
    bed_df = pd.read_table(bed_list[0], usecols=[0, 1, 2])
    bed_df.columns = ["chrom", "start", "end"]

    for bed_i in bed_list:
        logger.info(f"Load {bed_i} ...")
        sample_name = bed_i.parent.name
        df_i = pd.read_table(bed_i, usecols=[3])
        df_i.columns = [sample_name]
        if cov_cutoff is not None:
            if df_i[sample_name].quantile() < cov_cutoff:
                logger.info(f"跳过样本 {sample_name} (quantile < {cov_cutoff})")
                continue
        df_list.append(df_i)

    if not df_list:
        logger.error("所有样本都被过滤掉了，请检查 --cov-cutoff 参数")
        raise typer.Exit(code=1)

    df = pd.concat(df_list, axis=1)
    return bed_df, df


def get_stats_df(df: pd.DataFrame) -> pd.DataFrame:
    min_cov = df.min(axis=1)
    max_cov = df.max(axis=1)
    ave_cov = df.sum(axis=1) / df.shape[1]
    quantile_cov = df.quantile(0.5, axis=1)
    merged_df = pd.concat([min_cov, max_cov, ave_cov, quantile_cov], axis=1)
    merged_df.columns = ["min_cov", "max_cov", "mean_cov", "quantile_cov"]
    return merged_df


@app.command(help=MODULE_HELP)
def main(
    cds_cov_dir: Annotated[Path, typer.Argument(help="bamdst 输出目录（包含 */region.tsv.gz）")],
    out_file: Annotated[Path, typer.Argument(help="输出文件路径")],
    cov: Annotated[
        list[int] | None,
        typer.Option(help="覆盖度阈值列表"),
    ] = None,
    split_bed: Annotated[
        Path | None,
        typer.Option(help="split bed 文件路径（用于坐标转换）"),
    ] = None,
    cov_cutoff: Annotated[
        float | None,
        typer.Option(help="覆盖度过滤阈值（quantile 低于此值的样本将被跳过）"),
    ] = None,
    xlsx: Annotated[
        bool,
        typer.Option(help="同时输出 xlsx 格式文件"),
    ] = False,
) -> None:
    cov_list = cov if cov is not None else [1, 5, 10, 20, 30, 50, 100]
    bed_df, df_matrix = load_bed_files(cds_cov_dir, cov_cutoff=cov_cutoff)
    stats_df = get_stats_df(df_matrix)
    cov_df_list: list[pd.Series] = []

    for cov_i in cov_list:
        cov_i_df_matrix = df_matrix >= cov_i
        logger.info(f"Processing coverage {cov_i}x ...")
        cover_df = cov_i_df_matrix.sum(axis=1)
        cover_ratio_df = cover_df / cov_i_df_matrix.shape[1]
        cover_ratio_df.name = f"coverage_{cov_i}x"
        cov_df_list.append(cover_ratio_df)

    if split_bed is not None:
        bed_df = merge_chr(bed_df, split_bed)

    bed_df = bed_df[["chrom", "end"]].copy()
    bed_df = bed_df.rename(columns={"end": "pos"})

    result_df = pd.concat([bed_df, stats_df, *cov_df_list], axis=1)

    result_df.to_csv(out_file, index=False, float_format="%.3f", sep="\t")
    logger.success(f"TSV output written to: {out_file}")

    if xlsx:
        xlsx_path = out_file.with_suffix(".xlsx")
        result_df.to_excel(xlsx_path, index=False, float_format="%.3f")
        logger.success(f"XLSX output written to: {xlsx_path}")


if __name__ == "__main__":
    app()
