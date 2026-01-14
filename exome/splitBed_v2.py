"""
BED/FAI file splitting utility for exome analysis workflows.

This module provides functionality to split BED files or reference genome
index (.fai) files into multiple smaller files for parallel processing.
"""

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, cast

import polars as pl
import typer
from loguru import logger

OUT_COLUMNS = ["chrom", "start", "end"]

app = typer.Typer(
    help="Split BED files or FASTA index files for parallel processing.",
    no_args_is_help=True,
)


@dataclass
class BedRegion:
    """Represents a genomic region from a BED file."""

    chrom: str
    start: int
    end: int

    @property
    def length(self) -> int:
        """Return the length of the region."""
        return self.end - self.start


def calculate_padding(value: int) -> int:
    """Calculate the number of digits needed to represent a value."""
    return math.ceil(math.log10(value + 1)) if value > 0 else 1


def get_optimal_split_length(genome_length: int, split_number: int) -> int:
    """
    Calculate an optimal split length for genome splitting.

    The result is rounded to a nice number for easier processing tracking.

    Args:
        genome_length: Total genome length in base pairs.
        split_number: Number of desired splits.

    Returns:
        Optimal split length rounded to a meaningful value.
    """
    raw_length = genome_length // split_number
    magnitude = 10 ** int(math.log10(raw_length))
    multiplier = raw_length // magnitude
    return multiplier * magnitude


def build_output_filename(
    out_dir: Path,
    prefix_idx: str,
    start_region: BedRegion,
    end_region: BedRegion,
    pad_num: int,
) -> Path:
    """Generate output filename based on genomic coordinates."""
    start_pos = f"{start_region.chrom}_{start_region.start}"
    end_pos = str(end_region.end)

    if start_region.chrom != end_region.chrom:
        end_pos = f"{end_region.chrom}_{end_region.end}"

    return out_dir / f"{prefix_idx}_{start_pos}_{end_pos}.bed"


def save_bed_regions(
    regions: list[BedRegion],
    out_file: Path,
) -> None:
    """Save a list of BED regions to a file."""
    if not regions:
        return

    df = pl.DataFrame(
        [{"chrom": r.chrom, "start": r.start, "end": r.end} for r in regions],
        schema=OUT_COLUMNS,
    )
    df.write_csv(out_file, separator="\t", include_header=False)


def split_bed_file(
    bed_file: Path,
    out_dir: Path,
    split_number: int,
) -> None:
    """
    Split a BED file into multiple smaller files.

    The file is split based on total region length, aiming to create
    roughly equal-sized output files.

    Args:
        bed_file: Path to the input BED file.
        out_dir: Directory to write split files.
        split_number: Number of files to split into.
    """
    logger.info(f"Reading BED file: {bed_file}")

    bed_df = pl.read_csv(
        bed_file,
        separator="\t",
        has_header=False,
        new_columns=OUT_COLUMNS,
        schema_overrides={"start": pl.Int64, "end": pl.Int64},
    )

    bed_df = bed_df.with_columns((pl.col("end") - pl.col("start")).alias("region_length"))

    total_length = bed_df["region_length"].sum()
    if total_length is None or total_length == 0:
        logger.warning("Empty BED file")
        return

    length_per_file = total_length // split_number
    max_position = bed_df["end"].max()
    pad_num = calculate_padding(cast(int, max_position)) if max_position is not None else 1
    prefix_pad_num = calculate_padding(split_number) + 1

    output_dir = out_dir / bed_file.stem
    if output_dir.exists():
        logger.error(f"Output directory already exists: {output_dir}")
        raise typer.Exit(1)
    output_dir.mkdir(parents=True)

    logger.info(f"Splitting into ~{split_number} files, ~{length_per_file:,} bp per file")

    current_regions: list[BedRegion] = []
    current_size = 0
    current_idx = 0

    for row in bed_df.iter_rows(named=True):
        region = BedRegion(
            chrom=row["chrom"],
            start=row["start"],
            end=row["end"],
        )

        if current_size > length_per_file and current_regions:
            current_idx += 1
            prefix = str(current_idx).zfill(prefix_pad_num)
            out_file = build_output_filename(output_dir, prefix, current_regions[0], current_regions[-1], pad_num)
            save_bed_regions(current_regions, out_file)
            current_regions = []
            current_size = 0

        current_size += region.length
        current_regions.append(region)

    if current_regions:
        current_idx += 1
        prefix = str(current_idx).zfill(prefix_pad_num)
        out_file = build_output_filename(output_dir, prefix, current_regions[0], current_regions[-1], pad_num)
        save_bed_regions(current_regions, out_file)

    logger.success(f"Created {current_idx} split files in {output_dir}")


def split_fai_file(
    fai_file: Path,
    out_dir: Path,
    split_number: int,
) -> None:
    """
    Split a FASTA index (.fai) file into genomic chunks.

    Creates BED-like files representing chromosome segments for
    parallel processing of whole genomes.

    Args:
        fai_file: Path to the FASTA index file.
        out_dir: Directory to write split files.
        split_number: Number of files to split into.
    """
    logger.info(f"Reading FAI file: {fai_file}")

    fai_df = pl.read_csv(
        fai_file,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "chrom_length"],
        schema_overrides={"chrom_length": pl.Int64},
    )

    genome_length = fai_df["chrom_length"].sum()
    if genome_length is None or genome_length == 0:
        logger.warning("Empty FAI file")
        return

    split_length = get_optimal_split_length(cast(int, genome_length), split_number)
    max_chrom_length = fai_df["chrom_length"].max()
    pad_num = calculate_padding(cast(int, max_chrom_length)) if max_chrom_length is not None else 1
    prefix_pad_num = calculate_padding(split_number) + 1

    output_dir = out_dir / "genome"
    if output_dir.exists():
        logger.error(f"Output directory already exists: {output_dir}")
        raise typer.Exit(1)
    output_dir.mkdir(parents=True)

    logger.info(f"Splitting genome (~{genome_length:,} bp) into ~{split_number} files, ~{split_length:,} bp per file")

    step = split_length // 10
    current_regions: list[BedRegion] = []
    current_length = 0
    current_idx = 0

    for row in fai_df.iter_rows(named=True):
        chrom = row["chrom"]
        chrom_length = row["chrom_length"]

        for start in range(0, chrom_length, step):
            if current_length >= split_length and current_regions:
                current_idx += 1
                prefix = str(current_idx).zfill(prefix_pad_num)
                out_file = build_output_filename(
                    output_dir,
                    prefix,
                    current_regions[0],
                    current_regions[-1],
                    pad_num,
                )
                save_bed_regions(current_regions, out_file)
                current_regions = []
                current_length = 0

            end = min(start + step, chrom_length)
            region_length = end - start
            current_length += region_length

            # Merge with previous region if same chromosome
            if current_regions and current_regions[-1].chrom == chrom:
                last = current_regions[-1]
                current_regions[-1] = BedRegion(chrom=last.chrom, start=last.start, end=end)
            else:
                current_regions.append(BedRegion(chrom=chrom, start=start, end=end))

    if current_regions:
        current_idx += 1
        prefix = str(current_idx).zfill(prefix_pad_num)
        out_file = build_output_filename(output_dir, prefix, current_regions[0], current_regions[-1], pad_num)
        save_bed_regions(current_regions, out_file)

    logger.success(f"Created {current_idx} split files in {output_dir}")


@app.command()
def main(
    input_file: Annotated[
        Path,
        typer.Argument(help="Input BED or FAI file to split"),
    ],
    output_dir: Annotated[
        Path,
        typer.Argument(help="Output directory for split files"),
    ],
    split_number: Annotated[
        int,
        typer.Option(
            "--split-number",
            "-n",
            help="Number of files to split into",
            min=1,
        ),
    ] = 400,
    is_bed: Annotated[
        bool,
        typer.Option(
            "--bed/--fai",
            help="Treat input as BED file (default) or FAI file",
        ),
    ] = True,
) -> None:
    """
    Split a BED file or FASTA index file into multiple smaller files.

    For BED files: Splits regions based on total region length.

    For FAI files: Creates genomic coordinate chunks suitable for
    parallel processing of whole genomes.
    """
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        raise typer.Exit(1)

    try:
        if is_bed:
            split_bed_file(input_file, output_dir, split_number)
        else:
            split_fai_file(input_file, output_dir, split_number)
    except Exception as e:
        logger.exception(f"Error during splitting: {e}")
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
