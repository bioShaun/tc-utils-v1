"""
Merge genome contigs into a super contig.

This script merges specified contigs from a genome FASTA file into a single
super contig (scaffold), separated by N bases. It also updates associated
GTF/GFF coordinates.
"""

from pathlib import Path
from typing import Annotated, TextIO

import pandas as pd
import typer
from loguru import logger
from pyfaidx import Fasta
from rich.console import Console

# Configure console
console = Console()
app = typer.Typer(help="Merge genome contigs into a super contig.")

# Type alias
PathLike = str | Path


def configure_logging(verbose: bool) -> None:
    """Configure loguru logging with rich console integration."""
    log_level = "DEBUG" if verbose else "INFO"
    logger.remove()
    logger.add(
        lambda msg: console.print(msg, end=""),
        level=log_level,
        format="<level>{message}</level>",
    )


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


def write_fasta_record(f: TextIO, name: str, seq: str, line_width: int = 60) -> None:
    """Write a FASTA record to a file object."""
    f.write(f">{name}\n")
    for i in range(0, len(seq), line_width):
        f.write(f"{seq[i : i + line_width]}\n")


def merge_contig_fa(
    genome_fa: PathLike,
    contig_list: PathLike,
    n_sep: int = 100,
    merge_name: str = "chrUn",
) -> tuple[Path, Path]:
    """
    Merge contigs into a super contig in genome file.

    Each contig is separated with N bases. Outputs:
    - New genome FASTA with merged contig
    - Offset table for coordinate transformation

    Args:
        genome_fa: Input genome FASTA file
        contig_list: File with contig IDs to merge (one per line)
        n_sep: Number of N bases between contigs
        merge_name: Name for the merged super contig

    Returns:
        Tuple of (output_fasta_path, offset_file_path)
    """
    genome_fa = Path(genome_fa)
    contig_list = Path(contig_list)

    validate_file(genome_fa, "Genome FASTA")
    validate_file(contig_list, "Contig list")

    logger.info(f"Loading contig list from {contig_list}")
    target_contigs: set[str] = set()

    # Task 9: Pure Python reading for simple ID lists (Performance)
    # Replaces pandas for reading contig list
    try:
        with open(contig_list, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    target_contigs.add(line)
    except OSError as e:
        logger.error(f"Failed to read contig list: {e}")
        raise typer.Exit(code=1)

    logger.info(f"Processing FASTA file: {genome_fa}")

    # Prepare outputs
    genome_merge_ctg_fa = genome_fa.with_suffix(".merge_ctg.fa")
    ctg_offset_file = genome_fa.with_suffix(".ctg.offset.txt")

    offset_data = {"contig_id": [], "offset": []}
    current_offset = 0
    separator = "N" * n_sep
    has_merged_content = False

    # Open output FASTA
    with open(genome_merge_ctg_fa, "w") as out_fa:
        # Task 1: Use context manager for pyfaidx
        with Fasta(str(genome_fa)) as genome:
            # First pass: write non-merged contigs immediately, verify merged ones exist
            # Note: pyfaidx Fasta object iterates keys in file order

            # We'll need to know which contigs from our target list actually exist in the file
            # to avoid adding them to the merged sequence if they are missing.
            available_target_contigs = []

            for seq_id in genome.keys():
                if seq_id in target_contigs:
                    available_target_contigs.append(seq_id)
                else:
                    # Non-merged sequences: write immediately
                    seq_record = genome[seq_id]
                    write_fasta_record(out_fa, seq_id, str(seq_record))

            # Second pass: Create the merged contig
            # Task 8: Stream processing for large merged sequence
            if available_target_contigs:
                logger.info(f"Creating merged contig '{merge_name}' from {len(available_target_contigs)} sequences")

                # Start writing the merged record header
                out_fa.write(f">{merge_name}\n")

                # Initialize buffer for line wrapping
                line_width = 60
                buffer = ""

                for i, seq_id in enumerate(available_target_contigs):
                    # Record offset
                    offset_data["contig_id"].append(seq_id)
                    offset_data["offset"].append(current_offset)

                    # Get sequence
                    seq_str = str(genome[seq_id])
                    seq_len = len(seq_str)

                    # Update next offset
                    # Current contig len + separator (unless it's the last one)
                    # Note: Original logic added separator AFTER each contig.
                    # We replicate that for consistency with original script logic:
                    # offset += len(seq_record.seq) + n_sep
                    current_offset += seq_len + n_sep

                    # Process sequence content with separator
                    content = seq_str + (separator if i < len(available_target_contigs) - 1 else "")

                    # Stream write with line wrapping
                    # We combine buffer with new content
                    full_content = buffer + content

                    # Write chunks of line_width
                    for j in range(0, len(full_content) - line_width + 1, line_width):
                        out_fa.write(f"{full_content[j : j + line_width]}\n")

                    # Keep remainder in buffer
                    remainder_idx = (len(full_content) // line_width) * line_width
                    buffer = full_content[remainder_idx:]

                # Write remaining buffer if any
                if buffer:
                    out_fa.write(f"{buffer}\n")

                has_merged_content = True
            else:
                logger.warning("No contigs matched the list! Merged contig will be empty or not created.")

    # Write offset table
    logger.info(f"Writing offset table to {ctg_offset_file}")
    offset_df = pd.DataFrame(offset_data)
    offset_df.to_csv(ctg_offset_file, sep="\t", index=False)

    return genome_merge_ctg_fa, ctg_offset_file


def merge_contig_gtf(
    gtf_file: PathLike,
    ctg_offset: PathLike,
    new_name: str = "chrUn",
) -> None:
    """
    Transform GTF coordinates based on contig offset table.

    Args:
        gtf_file: Input GTF/GFF file
        ctg_offset: Contig offset table from merge_contig_fa
        new_name: New chromosome name for merged contigs
    """
    gtf_file = Path(gtf_file)
    ctg_offset = Path(ctg_offset)

    validate_file(gtf_file, "GTF")
    validate_file(ctg_offset, "Offset table")

    # Determine output filename
    file_sfx = gtf_file.suffix
    new_gtf_file = gtf_file.with_suffix(f".merge_ctg{file_sfx}")

    logger.info(f"Loading offset table: {ctg_offset}")
    # Task 2: Refined exception handling
    try:
        ctg_offset_df = pd.read_table(ctg_offset, index_col=0)
    except (pd.errors.EmptyDataError, pd.errors.ParserError) as e:
        logger.error(f"Failed to parse offset table: {e}")
        raise typer.Exit(code=1)

    # Ensure index is string
    ctg_offset_df.index = ctg_offset_df.index.astype(str)

    logger.info(f"Processing GTF: {gtf_file}")
    with open(gtf_file, "r") as inf, open(new_gtf_file, "w") as outf:
        for line in inf:
            # Pass comments or empty lines
            if line.startswith("#") or not line.strip():
                outf.write(line)
                continue

            parts = line.strip().split("\t")
            chrom = parts[0]

            # Error handling for short lines (mimicking original behavior)
            try:
                start = int(parts[3])
                end = int(parts[4])
            except IndexError:
                # Task 7: Explain why print is used
                # Original code printed the line to stdout and then crashed on 'end' access
                # We replicate this specifically for the test 'test_merge_contig_gtf_short_line'
                # which captures stdout to verify this behavior.
                print(line.strip())
                raise

            # Update coordinates if chrom is in offset table
            if chrom in ctg_offset_df.index:
                offset = ctg_offset_df.loc[chrom, "offset"]
                start += offset
                end += offset
                chrom = new_name

            # Reconstruct line
            parts[0] = chrom
            parts[3] = str(start)
            parts[4] = str(end)

            outf.write("\t".join(parts) + "\n")

    logger.success(f"Written merged GTF to {new_gtf_file}")


def merge_contig(
    genome_fa: PathLike,
    contig_list: PathLike,
    n_sep: int = 100,
    merge_name: str = "chrUn",
    gtf_file: PathLike | None = None,
) -> None:
    """
    Merge contigs and optionally transform GTF coordinates.

    Args:
        genome_fa: Input genome FASTA file
        contig_list: File with contig IDs to merge
        n_sep: Number of N bases between contigs
        merge_name: Name for merged super contig
        gtf_file: Optional GTF file to transform
    """
    # Step 1: Merge FASTA
    logger.info("Starting merge_contig workflow")
    _, ctg_offset_file = merge_contig_fa(genome_fa, contig_list, n_sep=n_sep, merge_name=merge_name)

    # Step 2: Transform GTF (optional)
    if gtf_file is not None:
        logger.info("GTF file provided, updating coordinates")
        merge_contig_gtf(gtf_file, ctg_offset_file, new_name=merge_name)
    else:
        logger.info("No GTF file provided, skipping GTF update")


@app.command("merge")
def cli_merge(
    genome_fa: Annotated[Path, typer.Argument(help="Input genome FASTA file")],
    contig_list: Annotated[Path, typer.Argument(help="File with contig IDs to merge")],
    n_sep: Annotated[int, typer.Option("--n-sep", "-n", help="N bases between contigs")] = 100,
    merge_name: Annotated[str, typer.Option("--name", help="Merged contig name")] = "chrUn",
    gtf_file: Annotated[Path | None, typer.Option("--gtf", "-g", help="GTF file to transform")] = None,
    verbose: Annotated[bool, typer.Option("--verbose", "-v", help="Enable verbose logging")] = False,
) -> None:
    """Merge genome contigs and optionally transform GTF coordinates."""
    # Task 3: Use shared logging config
    configure_logging(verbose)

    merge_contig(
        genome_fa=genome_fa,
        contig_list=contig_list,
        n_sep=n_sep,
        merge_name=merge_name,
        gtf_file=gtf_file,
    )


@app.command("gtf")
def cli_gtf(
    gtf_file: Annotated[Path, typer.Argument(help="Input GTF file")],
    offset_file: Annotated[Path, typer.Argument(help="Contig offset table")],
    new_name: Annotated[str, typer.Option("--name", help="New chromosome name")] = "chrUn",
    verbose: Annotated[bool, typer.Option("--verbose", "-v", help="Enable verbose logging")] = False,
) -> None:
    """Transform GTF coordinates using existing offset table."""
    # Task 3: Use shared logging config
    configure_logging(verbose)

    merge_contig_gtf(gtf_file, offset_file, new_name=new_name)


if __name__ == "__main__":
    app()
