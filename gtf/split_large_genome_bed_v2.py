"""
Split BED file coordinates based on genome split regions.

This script takes a BED file and a split regions file, then transforms the
coordinates from the original genome space to the split genome space.
"""

from pathlib import Path
from typing import Annotated

import polars as pl
import typer
from loguru import logger
from rich.console import Console

console = Console()
app = typer.Typer(help="Split BED file coordinates based on genome split regions.")


def validate_file(path: Path, name: str) -> None:
    """Validate that a file exists and is a file."""
    if not path.exists():
        logger.error(f"{name} file not found: {path}")
        console.print(f"[red]Error:[/red] {name} file not found: {path}")
        raise typer.Exit(code=1)

    if not path.is_file():
        logger.error(f"Path is not a file: {path}")
        console.print(f"[red]Error:[/red] Path is not a file: {path}")
        raise typer.Exit(code=1)


def transform_coordinates(
    bed_df: pl.DataFrame,
    split_df: pl.DataFrame,
) -> pl.DataFrame:
    """
    Transform BED coordinates from original to split genome space.

    Uses an efficient join strategy to map coordinates.
    """
    logger.info("Mapping BED regions to split regions")

    # Sort for join_asof - must be sorted by (chrom, start/split_start)
    bed_df = bed_df.sort(["chrom", "start"])
    split_df = split_df.sort(["chrom", "split_start"])

    # Explicitly mark columns as sorted to avoid Polars warning
    # This tells Polars the data is sorted within each 'chrom' group
    bed_df = bed_df.set_sorted("start")
    split_df = split_df.set_sorted("split_start")

    # join_asof finds the split where split_start <= start
    # Using backward strategy to match each region to its containing split
    result = bed_df.join_asof(
        split_df,
        left_on="start",
        right_on="split_start",
        by="chrom",
        strategy="backward",
    )

    # Filter regions that actually fall within the split boundaries
    # start must be >= split_start (guaranteed by join_asof backward)
    # and start must be < split_end
    filtered_df = result.filter(pl.col("split_start").is_not_null() & (pl.col("start") < pl.col("split_end")))

    initial_count = bed_df.height
    filtered_count = filtered_df.height

    if filtered_count < initial_count:
        dropped = initial_count - filtered_count
        logger.warning(f"Dropped {dropped} regions that fall outside split boundaries")

    if filtered_df.is_empty():
        logger.warning("No regions found within split boundaries")
        return filtered_df

    logger.info("Calculating new coordinates")
    # New coordinates are relative to the split start
    # Note: If a region spans multiple splits, this logic assigns it to the split where it starts.
    # The end coordinate might exceed the split length.
    filtered_df = filtered_df.with_columns(
        (pl.col("start") - pl.col("split_start")).alias("new_start"),
        (pl.col("end") - pl.col("split_start")).alias("new_end"),
    )

    return filtered_df


@app.command()
def main(
    bed: Annotated[Path, typer.Argument(help="Input BED file path")],
    split_bed: Annotated[Path, typer.Argument(help="Split regions BED file path")],
    output: Annotated[
        Path | None,
        typer.Option("--output", "-o", help="Output BED file path (default: <input>.split.bed)"),
    ] = None,
    verbose: Annotated[bool, typer.Option("--verbose", "-v", help="Enable verbose logging")] = False,
) -> None:
    """
    Split BED file coordinates based on genome split regions.

    This tool transforms coordinates from the original genome space to split
    genome space. It reads a BED file and a split regions file, then outputs
    a new BED file with transformed coordinates.

    The split regions file should have 4 columns:
    - chrom: chromosome name
    - split_start: start position of split region
    - split_end: end position of split region
    - new_chrom: new chromosome name for this region
    """
    # Configure logging
    log_level = "DEBUG" if verbose else "INFO"
    logger.remove()
    logger.add(
        lambda msg: console.print(msg, end=""),
        level=log_level,
        format="<level>{message}</level>",
    )

    # Validate inputs
    validate_file(bed, "BED")
    validate_file(split_bed, "Split BED")

    # Load data
    logger.info(f"Loading BED file: {bed}")
    # Read all columns, but we only know the first 3
    # We use has_header=False and then rename first columns
    bed_df = pl.read_csv(bed, separator="\t", has_header=False)
    # Rename first 3 columns, leave others as column_4, column_5 etc.
    new_names = {
        bed_df.columns[0]: "chrom",
        bed_df.columns[1]: "start",
        bed_df.columns[2]: "end",
    }
    bed_df = bed_df.rename(new_names)

    logger.info(f"Loading split BED file: {split_bed}")
    split_df = pl.read_csv(
        split_bed,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "split_start", "split_end", "new_chrom"],
    )

    # Transform coordinates
    result_df = transform_coordinates(bed_df, split_df)

    # Determine output path
    output_path = output if output else bed.with_suffix(".split.bed")

    # Save output
    if not result_df.is_empty():
        logger.info(f"Writing output to: {output_path}")

        # Construct output columns: new_chrom, new_start, new_end, then the rest
        other_cols = [c for c in result_df.columns if c.startswith("column_")]
        output_cols = ["new_chrom", "new_start", "new_end"] + other_cols

        result_df.select(output_cols).write_csv(output_path, separator="\t", include_header=False)
        logger.success(f"Successfully wrote {len(result_df)} regions to {output_path}")
    else:
        logger.warning("No output generated as no regions matched")


if __name__ == "__main__":
    app()
