"""
Extract representative (longest) transcript per gene from FASTA.

This script processes a FASTA file (typically from gffread) and extracts the
longest transcript sequence for each gene. It resolves gene IDs either from
FASTA headers or an external GTF/mapping file.
"""

import re
from pathlib import Path
from typing import Annotated

import polars as pl
import typer
from Bio import Seq, SeqIO
from loguru import logger
from rich.console import Console

# Configure console
console = Console()
app = typer.Typer(help="Extract representative (longest) transcript per gene from FASTA.")

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


def parse_gtf_attributes(attributes: str, key: str) -> str | None:
    """Extract value for a specific key from GTF attributes string."""
    # Regex to match key "value"; or key "value"
    match = re.search(f'{key} "([^"]+)"', attributes)
    if match:
        return match.group(1)
    return None


def load_transcript_gene_map(
    gtf: Path | None = None, 
    gene_tr_map: Path | None = None
) -> dict[str, str] | None:
    """
    Load transcript to gene mapping from GTF or map file.
    
    Returns:
        Dictionary mapping transcript_id -> gene_id, or None if no inputs provided.
    """
    if gtf:
        validate_file(gtf, "GTF")
        logger.info(f"Loading GTF file: {gtf}")
        try:
            # Read GTF with polars
            # GTF is tab-separated, no header, comments start with #
            df = pl.read_csv(
                gtf,
                separator="\t",
                has_header=False,
                comment_prefix="#",
                new_columns=[
                    "seqname", "source", "feature", "start", "end", 
                    "score", "strand", "frame", "attribute"
                ],
                schema_overrides={
                    "start": pl.Int64, 
                    "end": pl.Int64,
                    "score": pl.String, # score can be "."
                    "frame": pl.String  # frame can be "."
                }
            )
            
            # Filter for exons and select necessary columns
            df = df.filter(pl.col("feature") == "CDS").select("attribute")
            
            # Extract transcript_id and gene_id
            # We use map_elements for regex extraction as it's robust for complex GTF attributes
            # Note: For very large GTFs, this might be slower than pure regex on file line-by-line,
            # but polars provides a nice API. Given "exon" filter, data size is reduced.
            
            mapping = {}
            # Iterating rows for extraction (polars str.extract works but can be tricky with variable spacing)
            # Let's try a regex approach with polars if possible, or fallback to python iteration
            
            # Improved regex approach using polars string expressions
            df_extracted = df.with_columns([
                pl.col("attribute").str.extract(r'transcript_id "([^"]+)"', 1).alias("transcript_id"),
                pl.col("attribute").str.extract(r'gene_id "([^"]+)"', 1).alias("gene_id")
            ]).filter(
                pl.col("transcript_id").is_not_null() & pl.col("gene_id").is_not_null()
            ).select(["transcript_id", "gene_id"]).unique()
            
            # Convert to dictionary
            return dict(zip(df_extracted["transcript_id"], df_extracted["gene_id"]))

        except Exception as e:
            logger.error(f"Failed to parse GTF with polars: {e}")
            raise typer.Exit(code=1)

    elif gene_tr_map:
        validate_file(gene_tr_map, "Gene-Transcript Map")
        logger.info(f"Loading map file: {gene_tr_map}")
        try:
            # Expecting format: transcript_id <tab> gene_id OR gene_id <tab> transcript_id?
            # Original code: index_col=1, names=["gene_id"] -> implies col 0 is index (transcript_id?), col 1 is gene_id?
            # Original: pd.read_csv(gene_tr_map, header=None, index_col=1, names=["gene_id"], sep="\t")
            # If names=["gene_id"], and we set it, pandas usually takes the remaining cols.
            # Let's assume standard format: transcript_id\tgene_id
            
            # Re-reading original logic:
            # pd.read_csv(gene_tr_map, header=None, index_col=1, names=["gene_id"], sep="\t")
            # If input is: T1\tG1
            # read_csv(..., names=["gene_id"]) -> usually names applies to ALL columns if header=None
            # This looks suspicious. Let's assume input is 2 columns.
            # If index_col=1, then column 1 is index.
            # If the file has 2 columns: Col0, Col1. Index is Col1.
            # So map is: Col1 (Gene) -> Col0 (Transcript)?
            # But the usage later is: `if seq_record.id in tr_gtf_df.index: gene_id = ...`
            # So index MUST be transcript_id.
            # So Col1 must be transcript_id?
            
            # Let's be robust: Read 2 columns, user can specify if needed, or assume col 0 is transcript, col 1 is gene
            # Standard gffread output or simple maps are usually: gene_id transcript_id OR transcript_id gene_id
            
            # Let's replicate original pandas behavior exactly to be safe, then optimize.
            # Original: index_col=1. So 2nd column is index (transcript_id).
            # So input file structure expected: gene_id <tab> transcript_id
            
            df = pl.read_csv(
                gene_tr_map, 
                separator="\t", 
                has_header=False,
                new_columns=["gene_id", "transcript_id"] # Assuming this order based on index_col=1
            )
            
            return dict(zip(df["transcript_id"], df["gene_id"]))
            
        except Exception as e:
            logger.error(f"Failed to parse map file: {e}")
            raise typer.Exit(code=1)

    return None


def extract_gene_id(
    description: str, 
    seq_id: str, 
    tr_map: dict[str, str] | None
) -> str | None:
    """Extract gene_id from sequence description or mapping."""
    if tr_map:
        if seq_id in tr_map:
            return tr_map[seq_id]
        else:
            return None # Will trigger skip in main loop

    # Fallback to regex on description
    # Pattern: gene=GENE_ID
    match = re.search(r"gene=(\S+)", description)
    if match:
        return match.group(1)
    
    logger.error(f"Wrong fasta format. [>transcript_id gene=gene_id] not found in: {description}")
    raise typer.Exit(code=1)


def gene_presentative_fa(
    gffread_fa: PathLike,
    gtf: PathLike | None = None,
    gene_tr_map: PathLike | None = None,
) -> None:
    """
    Extract representative (longest) transcript per gene from FASTA.

    Args:
        gffread_fa: Input FASTA file from gffread
        gtf: GTF file for transcript-gene mapping
        gene_tr_map: Alternative TSV file for transcript-gene mapping

    Raises:
        typer.Exit: If input validation fails or gene ID cannot be extracted
    """
    gffread_fa = Path(gffread_fa)
    if gtf: gtf = Path(gtf)
    if gene_tr_map: gene_tr_map = Path(gene_tr_map)

    validate_file(gffread_fa, "Input FASTA")
    
    # Load mapping
    tr_map = load_transcript_gene_map(gtf, gene_tr_map)
    
    logger.info(f"Processing FASTA: {gffread_fa}")
    
    # Dictionary to store longest transcript per gene
    # Format: gene_id -> (length, record)
    gene_best_tr: dict[str, tuple[int, SeqIO.SeqRecord]] = {}
    
    # Output file
    output_fa = gffread_fa.with_suffix(".gene.fa")
    
    try:
        # Biopython SeqIO is good here because we need to modify sequences and attributes
        # and memory usage is usually manageable for transcriptomes (vs whole genomes)
        for record in SeqIO.parse(gffread_fa, "fasta"):
            # Replace dot in seq (diamond mkindex fail compatibility)
            # Note: modifying record.seq in place
            if "." in record.seq:
                # Convert to mutable if needed, or create new Seq
                # Biopython Seq is immutable-ish in recent versions, strings are immutable
                record.seq = Seq.Seq(str(record.seq).replace(".", "*"))
            
            # Determine Gene ID
            gene_id = extract_gene_id(record.description, record.id, tr_map)
            
            if not gene_id:
                continue
                
            # Update record attributes
            seq_len = len(record.seq)
            record.id = gene_id
            record.description = "" # Clear description as requested
            
            # Check if this is the longest for this gene
            if gene_id in gene_best_tr:
                current_best_len, _ = gene_best_tr[gene_id]
                if seq_len <= current_best_len:
                    continue
            
            gene_best_tr[gene_id] = (seq_len, record)
            
    except Exception as e:
        logger.error(f"Error processing FASTA: {e}")
        raise typer.Exit(code=1)

    # Write output
    logger.info(f"Writing {len(gene_best_tr)} representative transcripts to {output_fa}")
    records_to_write = [data[1] for data in gene_best_tr.values()]
    
    if records_to_write:
        with open(output_fa, "w") as f:
            SeqIO.write(records_to_write, f, "fasta")
        logger.success("Done")
    else:
        logger.warning("No transcripts extracted!")


@app.command()
def main(
    gffread_fa: Annotated[Path, typer.Argument(help="Input FASTA file from gffread")],
    gtf: Annotated[Path | None, typer.Option("--gtf", "-g", help="GTF file for transcript-gene mapping")] = None,
    gene_tr_map: Annotated[Path | None, typer.Option("--map", "-m", help="Transcript-Gene mapping file")] = None,
    verbose: Annotated[bool, typer.Option("--verbose", "-v", help="Enable verbose logging")] = False,
) -> None:
    """Extract representative (longest) transcript per gene from FASTA."""
    configure_logging(verbose)
    gene_presentative_fa(gffread_fa, gtf, gene_tr_map)


if __name__ == "__main__":
    app()
