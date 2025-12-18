#!/usr/bin/env python3
"""
Merge two large tables (Table A and Table B) using Polars.
- Table A: Database table with id and variant information
- Table B: Annotation table with type, impact, gene, transcript, etc.

Join on: chrom, pos, refer, alt
Output: id (from A), type, impact, gene, transcript, exon_rank, cds_pos, protein_pos (from B)

Uses Polars lazy evaluation for memory-efficient processing of GB-level files.
"""

from pathlib import Path
from typing import Annotated

import polars as pl
import typer
from loguru import logger

app = typer.Typer(
    help="Merge Table A (with id) and Table B (with annotations) on variant keys"
)


def read_table_a_lazy(path: str, separator: str) -> pl.LazyFrame:
    """
    Read Table A lazily, selecting only necessary columns for the join.
    Columns needed: id, chrom, pos, refer, alt
    """
    return pl.scan_csv(
        path,
        separator=separator,
        has_header=True,
        low_memory=True,
    ).select(
        [
            pl.col("id"),
            pl.col("chrom"),
            pl.col("pos").cast(pl.Int64),
            pl.col("refer"),
            pl.col("alt"),
        ]
    )


def read_table_b_lazy(path: str, separator: str, has_header: bool) -> pl.LazyFrame:
    """
    Read Table B lazily.
    Expected columns: chrom, pos, refer, alt, type, impact, gene, transcript, exon_rank, cds_pos, protein_pos
    Note: exon_rank and protein_pos may be missing in some rows.
    """
    column_names = [
        "chrom",
        "pos",
        "refer",
        "alt",
        "type",
        "impact",
        "gene",
        "transcript",
        "exon_rank",
        "cds_pos",
        "protein_pos",
    ]

    if has_header:
        lf = pl.scan_csv(
            path,
            separator=separator,
            has_header=True,
            low_memory=True,
        )
    else:
        lf = pl.scan_csv(
            path,
            separator=separator,
            has_header=False,
            new_columns=column_names,
            low_memory=True,
        )

    # Select columns, handling missing ones with allow_missing=True
    return lf.select(
        pl.col("chrom"),
        pl.col("pos").cast(pl.Int64),
        pl.col("refer"),
        pl.col("alt"),
        pl.col("type"),
        pl.col("impact"),
        pl.col("gene"),
        pl.col("transcript"),
        pl.col("exon_rank").fill_null(""),
        pl.col("cds_pos"),
        pl.col("protein_pos").fill_null(""),
    )


def merge_tables(
    table_a: pl.LazyFrame,
    table_b: pl.LazyFrame,
) -> pl.LazyFrame:
    """
    Join Table A and Table B on chrom, pos, refer, alt.
    Output columns: id (from A), type, impact, gene, transcript, exon_rank, cds_pos, protein_pos (from B)
    """
    # Perform inner join on variant keys
    joined = table_a.join(
        table_b,
        on=["chrom", "pos", "refer", "alt"],
        how="inner",
    )

    # Select output columns
    return joined.select(
        [
            pl.col("id"),
            pl.col("type"),
            pl.col("impact"),
            pl.col("gene"),
            pl.col("transcript"),
            pl.col("exon_rank"),
            pl.col("cds_pos"),
            pl.col("protein_pos"),
        ]
    )


@app.command()
def main(
    table_a: Annotated[
        Path,
        typer.Option(
            "-a",
            "--table-a",
            help="Path to Table A (TSV with header, contains id column)",
        ),
    ],
    table_b: Annotated[
        Path,
        typer.Option("-b", "--table-b", help="Path to Table B (TSV, annotation table)"),
    ],
    output: Annotated[
        Path,
        typer.Option("-o", "--output", help="Output file path (TSV format)"),
    ],
    table_b_header: Annotated[
        bool,
        typer.Option("--table-b-header", help="Table B has header line"),
    ] = False,
    sep_a: Annotated[
        str,
        typer.Option("--sep-a", help="Separator for Table A", show_default="tab"),
    ] = "\t",
    sep_b: Annotated[
        str,
        typer.Option("--sep-b", help="Separator for Table B", show_default="tab"),
    ] = "\t",
):
    """Merge Table A and Table B on variant keys (chrom, pos, refer, alt)."""
    # Validate input files
    if not table_a.exists():
        logger.error(f"Table A not found: {table_a}")
        raise typer.Exit(1)
    if not table_b.exists():
        logger.error(f"Table B not found: {table_b}")
        raise typer.Exit(1)

    logger.info(f"Reading Table A: {table_a}")
    lf_a = read_table_a_lazy(str(table_a), sep_a)

    logger.info(f"Reading Table B: {table_b}")
    lf_b = read_table_b_lazy(str(table_b), sep_b, table_b_header)

    logger.info("Merging tables on chrom, pos, refer, alt...")
    result = merge_tables(lf_a, lf_b)

    # Use streaming for memory-efficient output
    logger.info(f"Writing output to: {output}")
    result.sink_csv(str(output), separator="\t")

    logger.success("Done!")


if __name__ == "__main__":
    app()
