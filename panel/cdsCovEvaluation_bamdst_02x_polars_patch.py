import logging
import sys
from functools import reduce
from pathlib import Path
from typing import List, Optional, Tuple, Union

import polars as pl
import typer
from loguru import logger


# 将 loguru 日志转发到标准 logging
class InterceptHandler(logging.Handler):
    def emit(self, record):
        logging.getLogger(record.name).handle(record)


logger.remove()
logger.add(InterceptHandler(), format="{message}", level="DEBUG")

DEFAULT_DEPTH_THRESHOLD = 0.2
DEFAULT_FLOAT_PRECISION = 3


class CoverageAnalysisError(Exception):
    pass


class DataValidationError(CoverageAnalysisError):
    pass


class FileProcessingError(CoverageAnalysisError):
    pass


def validate_input_path(path: Path, path_type: str = "directory") -> None:
    if not path.exists():
        raise DataValidationError(f"{path_type.capitalize()} does not exist: {path}")
    if path_type == "directory" and not path.is_dir():
        raise DataValidationError(f"Path is not a directory: {path}")
    elif path_type == "file" and not path.is_file():
        raise DataValidationError(f"Path is not a file: {path}")


def validate_dataframe(
    df: pl.DataFrame, name: str, required_columns: Optional[List[str]] = None
) -> None:
    if df.is_empty():
        logging.warning(f"{name} is empty")
        return
    if required_columns:
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise DataValidationError(
                f"{name} missing required columns: {missing_columns}"
            )


def merge_chr(df: pl.DataFrame, split_bed: Path) -> pl.DataFrame:
    try:
        validate_input_path(split_bed, "file")
        logging.info(f"Reading split bed file: {split_bed}")
        split_bed_df = pl.read_csv(
            split_bed,
            has_header=False,
            new_columns=["new_chrom", "offset", "offset_end", "chrom"],
            separator="\t",
        )
        validate_dataframe(
            split_bed_df, "split_bed_df", ["new_chrom", "offset", "offset_end", "chrom"]
        )
        merged_df = df.join(split_bed_df, on="chrom", how="inner")
        if merged_df.is_empty():
            logging.warning("No matching chromosomes found in split bed file")
            return df
        # 计算新坐标并重命名
        merged_df = merged_df.with_columns(
            (pl.col("start") + pl.col("offset")).alias("start"),
            (pl.col("end") + pl.col("offset")).alias("end"),
        )
        final_cols = ["new_chrom", "start", "end"] + [
            col
            for col in merged_df.columns
            if col not in ["new_chrom", "start", "end", "offset", "offset_end", "chrom"]
        ]
        merged_df = merged_df.select(final_cols).rename({"new_chrom": "chrom"})
        logging.info(
            f"Chromosome coordinates transformed for {merged_df.shape[0]} regions"
        )
        return merged_df
    except pl.exceptions.PolarsError as e:
        raise FileProcessingError(f"Failed to process split bed file: {e}")
    except DataValidationError:
        raise
    except Exception as e:
        raise FileProcessingError(f"Unexpected error in merge_chr: {e}")


def load_single_depth_file(
    depth_file: Path, sample_name: str
) -> Optional[pl.DataFrame]:
    try:
        logging.debug(f"Loading depth file: {depth_file}")
        df_i = pl.read_csv(depth_file, separator="	", has_header=True)
        if df_i.is_empty():
            logging.warning(f"Empty depth file: {depth_file}")
            return None

        depth_col_name = df_i.columns[2]

        df_with_coords = df_i.select(
            [
                pl.col("#Chr").alias("chrom"),
                pl.col("Pos").alias("end"),
                pl.col(depth_col_name),
            ]
        )

        mean_depth = df_with_coords[depth_col_name].mean()
        if mean_depth is None or mean_depth == 0:
            logging.warning(f"Zero or null mean depth in {sample_name}")
            depth_threshold = 0
        else:
            depth_threshold = mean_depth * DEFAULT_DEPTH_THRESHOLD

        logging.debug(
            f"Sample {sample_name}: record_num={df_i.shape[0]}, mean_depth={mean_depth:.2f if mean_depth else 0}, threshold={depth_threshold:.2f}"
        )

        df_i = df_with_coords.with_columns(
            [
                (pl.col("end") - 1).alias("start"),
                (pl.col(depth_col_name) >= depth_threshold).alias(sample_name),
            ]
        ).select(["chrom", "start", "end", sample_name])
        return df_i
    except pl.exceptions.NoDataError:
        logging.warning(f"Empty depth file: {depth_file}")
        return None
    except Exception as e:
        logging.error(f"Failed to load {depth_file}: {e}")
        return None


def load_bed_files(
    bed_dir: Path,
    sample_list: Optional[List[str]] = None,
) -> pl.DataFrame:
    try:
        validate_input_path(bed_dir, "directory")
        bed_list = sorted(list(bed_dir.glob("*/depth.tsv.gz")))
        if not bed_list:
            raise FileProcessingError(f"No depth.tsv files found in {bed_dir}")
        logging.info(f"Found {len(bed_list)} depth files")

        df_list = []
        processed_samples = []
        for bed_file in bed_list:
            sample_name = bed_file.parent.name
            if sample_list is not None and sample_name not in sample_list:
                logging.debug(f"Skipping sample {sample_name} (not in sample list)")
                continue
            logging.info(f"Processing sample: {sample_name}")
            df_i = load_single_depth_file(bed_file, sample_name)
            if df_i is not None:
                df_list.append(df_i)
                processed_samples.append(sample_name)
            else:
                logging.warning(f"Failed to process sample: {sample_name}")
        if not df_list:
            logging.warning("No valid samples processed, returning empty dataframe")
            return pl.DataFrame()

        logging.info(f"Merging data from {len(df_list)} samples: {processed_samples}")
        merged_df = reduce(
            lambda left, right: left.join(
                right, on=["chrom", "start", "end"], how="outer"
            ),
            df_list,
        )
        # For regions not present in a sample, coverage is False (0).
        merged_df = merged_df.fill_null(False)

        logging.info(f"Final merged dataframe shape: {merged_df.shape}")
        return merged_df

    except (FileProcessingError, DataValidationError):
        raise
    except Exception as e:
        raise FileProcessingError(f"Unexpected error in load_bed_files: {e}")


def calculate_coverage_ratio(df_matrix: pl.DataFrame) -> pl.Series:
    if df_matrix.is_empty():
        return pl.Series("coverage_0.2x", [])
    cover_sum = df_matrix.sum_horizontal()
    cover_ratio = cover_sum / df_matrix.width
    logging.debug(f"Coverage ratio calculated for {len(cover_ratio)} regions")
    return cover_ratio


def load_sample_list(sample_path: Path) -> List[str]:
    try:
        validate_input_path(sample_path, "file")
        # 如果文件为空，直接返回 []
        if sample_path.stat().st_size == 0:
            return []
        sample_df = pl.read_csv(sample_path, has_header=False, new_columns=["sample"])
        if sample_df.is_empty():
            return []
        sample_list = sample_df["sample"].to_list()
        logging.info(f"Loaded {len(sample_list)} samples from {sample_path}")
        return sample_list
    except Exception as e:
        raise FileProcessingError(f"Failed to load sample list from {sample_path}: {e}")


def write_output(
    df: pl.DataFrame, output_file: Path, float_precision: int = DEFAULT_FLOAT_PRECISION
) -> None:
    try:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        try:
            df.write_csv(output_file, separator="\t", float_precision=float_precision)
        except TypeError:
            df.write_csv(output_file, separator="\t")
        logging.info(f"Results written to: {output_file}")
        logging.info(f"Output contains {df.shape[0]} regions and {df.shape[1]} columns")
    except Exception as e:
        raise FileProcessingError(f"Failed to write output to {output_file}: {e}")


def main(
    cds_cov_dir: Path,
    out_file: Path,
    split_bed: Optional[Path] = typer.Option(
        None, help="Split bed file for coordinate transformation"
    ),
    sample_path: Optional[Path] = typer.Option(
        None, help="File containing sample list (one per line)"
    ),
    log_level: str = typer.Option(
        "INFO", help="Log level (DEBUG, INFO, WARNING, ERROR)"
    ),
) -> None:
    logging.basicConfig(
        level=log_level.upper(), format="%(asctime)s | %(levelname)s | %(message)s"
    )
    try:
        logging.info("Starting CDS coverage analysis")
        logging.info(f"Input directory: {cds_cov_dir}")
        logging.info(f"Output file: {out_file}")
        sample_list = None
        if sample_path is not None:
            logging.info(f"Loading sample list from: {sample_path}")
            sample_list = load_sample_list(sample_path)
        logging.info("Loading BED files and depth data")
        merged_df = load_bed_files(cds_cov_dir, sample_list=sample_list)

        if merged_df.is_empty() or merged_df.width <= 3:
            logging.warning("No data loaded, creating empty output file")
            if merged_df.is_empty():
                bed_df = pl.DataFrame(
                    schema={"chrom": pl.Utf8, "start": pl.Int64, "end": pl.Int64}
                )
            else:
                bed_df = merged_df.select(["chrom", "start", "end"])
            cover_ratio_df = bed_df.with_columns(pl.lit(0.0).alias("coverage_0.2x"))
        else:
            logging.info("Calculating coverage ratios")
            bed_df = merged_df.select(["chrom", "start", "end"])
            df_matrix = merged_df.drop(["chrom", "start", "end"])
            cover_ratio = calculate_coverage_ratio(df_matrix)
            cover_ratio_df = pl.concat(
                [bed_df, pl.DataFrame({"coverage_0.2x": cover_ratio})],
                how="horizontal",
            )
        if split_bed is not None:
            logging.info(f"Applying coordinate transformation using: {split_bed}")
            cover_ratio_df = merge_chr(cover_ratio_df, split_bed)
        write_output(cover_ratio_df, out_file)
        logging.info("Analysis completed successfully")
    except (CoverageAnalysisError, typer.Exit) as e:
        logging.error(f"Analysis failed: {e}")
        raise typer.Exit(1)
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        raise typer.Exit(1)


if __name__ == "__main__":
    typer.run(main)
