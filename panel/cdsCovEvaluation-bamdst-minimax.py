"""
Coverage Evaluation Script for BAM Depth Files

This script analyzes coverage statistics from bamdst output files,
supporting both region-based and site-based depth analysis.
"""

from functools import reduce
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import typer
from loguru import logger

# Configuration constants
DEFAULT_COVERAGE_THRESHOLDS = [1, 5, 10, 20, 30, 50, 100]
REGION_FILE_PATTERN = "*/region.tsv.gz"
DEPTH_FILE_PATTERN = "*/depth.tsv.gz"
BED_COLUMNS = ["chrom", "start", "end", "transcript_id"]


class FileLoader:
    """Handles loading and processing of coverage data files."""

    @staticmethod
    def _check_files_exist(file_list: List[Path]) -> None:
        """Verify that all required files exist."""
        if not file_list:
            raise FileNotFoundError("No files found matching the pattern")
        for file_path in file_list:
            if not file_path.exists():
                raise FileNotFoundError(f"File not found: {file_path}")

    @staticmethod
    def _validate_dataframe(df: pd.DataFrame, required_cols: List[str]) -> None:
        """Validate that DataFrame contains required columns."""
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")

    @staticmethod
    def load_region_files(
        bed_dir: Path, cov_cutoff: Optional[float] = None
    ) -> Tuple[pd.DataFrame, List[pd.DataFrame]]:
        """
        Load coverage data from region.tsv.gz files.

        Args:
            bed_dir: Directory containing sample subdirectories with region.tsv.gz files
            cov_cutoff: Optional coverage cutoff to filter samples

        Returns:
            Tuple of (bed_dataframe, list_of_sample_dataframes)

        Raises:
            FileNotFoundError: If no region files are found
            ValueError: If files don't have expected structure
        """
        bed_list = sorted(list(bed_dir.glob(REGION_FILE_PATTERN)))
        FileLoader._check_files_exist(bed_list)

        logger.info(f"Found {len(bed_list)} region files")

        # Load BED coordinates from first file
        logger.info(f"Load {bed_list[0]} ...")
        bed_df = pd.read_table(bed_list[0], usecols=[0, 1, 2])
        bed_df.columns = ["chrom", "start", "end"]

        df_list = []
        for bed_file in bed_list:
            logger.info(f"Load {bed_file} ...")
            sample_name = bed_file.parent.name

            # Load coverage depth column
            df_i = pd.read_table(bed_file, usecols=[3])
            df_i.columns = [sample_name]

            # Apply coverage cutoff filter
            if cov_cutoff is not None:
                if df_i[sample_name].quantile() < cov_cutoff:
                    logger.info(
                        f"Skipping {sample_name}: quantile {df_i[sample_name].quantile():.2f} "
                        f"below cutoff {cov_cutoff}"
                    )
                    continue

            df_list.append(df_i)

        if not df_list:
            raise ValueError("No samples passed the coverage cutoff filter")

        return bed_df, df_list

    @staticmethod
    def load_depth_files(
        bed_dir: Path, cov_cutoff: Optional[float] = None
    ) -> Tuple[pd.DataFrame, List[pd.DataFrame]]:
        """
        Load coverage data from depth.tsv.gz files.

        Args:
            bed_dir: Directory containing sample subdirectories with depth.tsv.gz files
            cov_cutoff: Optional coverage cutoff to filter samples

        Returns:
            Tuple of (bed_dataframe, list_of_sample_dataframes)

        Raises:
            FileNotFoundError: If no depth files are found
            ValueError: If files don't have expected structure
        """
        depth_list = sorted(list(bed_dir.glob(DEPTH_FILE_PATTERN)))
        FileLoader._check_files_exist(depth_list)

        logger.info(f"Found {len(depth_list)} depth files")

        # Load BED coordinates from first file
        logger.info(f"Load {depth_list[0]} ...")
        bed_df = pd.read_table(depth_list[0], usecols=["#Chr", "Pos"])

        # Validate required columns
        FileLoader._validate_dataframe(bed_df, ["#Chr", "Pos"])

        # Convert to BED format: start = Pos - 1, end = Pos
        bed_df = pd.DataFrame(
            {"chrom": bed_df["#Chr"], "start": bed_df["Pos"] - 1, "end": bed_df["Pos"]}
        )

        df_list = []
        for depth_file in depth_list:
            logger.info(f"Load {depth_file} ...")
            sample_name = depth_file.parent.name

            # Load "Cover depth" column
            df_i = pd.read_table(depth_file, usecols=["Cover depth"])
            df_i.columns = [sample_name]

            # Apply coverage cutoff filter
            if cov_cutoff is not None:
                if df_i[sample_name].quantile() < cov_cutoff:
                    logger.info(
                        f"Skipping {sample_name}: quantile {df_i[sample_name].quantile():.2f} "
                        f"below cutoff {cov_cutoff}"
                    )
                    continue

            df_list.append(df_i)

        if not df_list:
            raise ValueError("No samples passed the coverage cutoff filter")

        return bed_df, df_list


class DataMerger:
    """Handles merging and combining DataFrames."""

    @staticmethod
    def merge_sample_dataframes(df_list: List[pd.DataFrame]) -> pd.DataFrame:
        """
        Merge multiple sample DataFrames into a single matrix.

        Args:
            df_list: List of DataFrames, each containing coverage data for one sample

        Returns:
            Merged DataFrame with samples as columns
        """
        if not df_list:
            raise ValueError("No DataFrames to merge")

        if len(df_list) == 1:
            logger.warning("Only one sample found in the dataset")
            return df_list[0]

        logger.info(f"Merging {len(df_list)} sample DataFrames")
        merged_df = reduce(
            lambda x, y: pd.merge(x, y, left_index=True, right_index=True),
            df_list,
        )
        return merged_df


class CoverageAnalyzer:
    """Analyzes coverage statistics and generates metrics."""

    @staticmethod
    def calculate_summary_stats(df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate summary statistics for coverage data.

        Args:
            df: DataFrame with samples as columns and genomic positions as rows

        Returns:
            DataFrame with summary statistics columns
        """
        logger.info("Calculating summary statistics")

        stats = {
            "min_cov": df.min(axis=1),
            "max_cov": df.max(axis=1),
            "mean_cov": df.mean(axis=1),
            "median_cov": df.median(axis=1),
        }

        stats_df = pd.DataFrame(stats)
        return stats_df

    @staticmethod
    def calculate_coverage_ratios(
        df: pd.DataFrame, thresholds: List[int]
    ) -> List[pd.Series]:
        """
        Calculate coverage ratios for multiple threshold values.

        Args:
            df: DataFrame with samples as columns and genomic positions as rows
            thresholds: List of coverage thresholds to evaluate

        Returns:
            List of Series, each containing coverage ratios for a threshold
        """
        logger.info(f"Calculating coverage ratios for thresholds: {thresholds}")

        cov_ratios = []
        num_samples = df.shape[1]

        for threshold in thresholds:
            logger.debug(f"Processing threshold: {threshold}x")
            # Boolean matrix: True if coverage >= threshold
            meets_threshold = df >= threshold
            # Count samples meeting threshold for each position
            covered_samples = meets_threshold.sum(axis=1)
            # Calculate ratio
            coverage_ratio = covered_samples / num_samples
            coverage_ratio.name = f"coverage_{threshold}x"

            cov_ratios.append(coverage_ratio)

        return cov_ratios


class CoordinateTransformer:
    """Handles coordinate transformations for split genome scenarios."""

    @staticmethod
    def merge_chr_coordinates(df: pd.DataFrame, split_bed: Path) -> pd.DataFrame:
        """
        Transform coordinates based on split genome information.

        Args:
            df: DataFrame with genomic coordinates
            split_bed: Path to split BED file containing offset information

        Returns:
            DataFrame with transformed coordinates
        """
        logger.info(f"Merging chromosome coordinates using {split_bed}")

        split_bed_df = pd.read_csv(
            split_bed,
            header=None,
            names=["new_chrom", "offset", "offset_end", "chrom"],
            sep="\t",
        )

        merged_df = df.merge(split_bed_df)
        merged_df["new_start"] = merged_df["start"] + merged_df["offset"]
        merged_df["new_end"] = merged_df["end"] + merged_df["offset"]

        # Clean up columns
        columns_to_drop = ["chrom", "offset", "offset_end", "start", "end"]
        merged_df.drop(columns_to_drop, axis=1, inplace=True)

        # Rename columns
        merged_df.rename(
            columns={
                "new_chrom": "chrom",
                "new_start": "start",
                "new_end": "end",
            },
            inplace=True,
        )

        return merged_df


def main(
    cds_cov_dir: Path,
    out_file: Path,
    cov: List[int] = typer.Option(DEFAULT_COVERAGE_THRESHOLDS),
    split_bed: Optional[Path] = typer.Option(None),
    cov_cutoff: Optional[float] = typer.Option(None),
    use_site: bool = typer.Option(
        False, help="Use depth.tsv.gz for analysis (site-based) instead of region.tsv.gz (region-based)"
    ),
) -> None:
    """
    Main entry point for coverage evaluation.

    Args:
        cds_cov_dir: Directory containing sample subdirectories with coverage files
        out_file: Output TSV file path for results
        cov: List of coverage thresholds to evaluate
        split_bed: Optional path to split BED file for coordinate transformation
        cov_cutoff: Optional coverage cutoff to filter samples
        use_site: If True, use depth.tsv.gz files; otherwise use region.tsv.gz files
    """
    # Load data based on mode
    if use_site:
        logger.info("Using site-based depth analysis (depth.tsv.gz files)")
        bed_df, sample_df_list = FileLoader.load_depth_files(
            cds_cov_dir, cov_cutoff=cov_cutoff
        )
    else:
        logger.info("Using region-based depth analysis (region.tsv.gz files)")
        bed_df, sample_df_list = FileLoader.load_region_files(
            cds_cov_dir, cov_cutoff=cov_cutoff
        )

    # Merge sample data
    coverage_matrix = DataMerger.merge_sample_dataframes(sample_df_list)

    # Calculate statistics
    stats_df = CoverageAnalyzer.calculate_summary_stats(coverage_matrix)
    coverage_ratios = CoverageAnalyzer.calculate_coverage_ratios(coverage_matrix, cov)

    # Apply coordinate transformation if needed
    if split_bed is not None:
        logger.info("Applying chromosome coordinate transformation")
        bed_df = CoordinateTransformer.merge_chr_coordinates(bed_df, split_bed)

    # Combine all results
    logger.info("Combining results")
    final_df = pd.concat([bed_df, stats_df, *coverage_ratios], axis=1)

    # Write output
    logger.info(f"Writing results to {out_file}")
    out_file.parent.mkdir(parents=True, exist_ok=True)
    final_df.to_csv(out_file, index=False, float_format="%.3f", sep="\t")

    # Print summary
    num_samples = len(sample_df_list)
    num_positions = len(final_df)
    logger.success(
        f"Analysis complete: {num_samples} samples, {num_positions} positions, "
        f"output saved to {out_file}"
    )


if __name__ == "__main__":
    typer.run(main)
