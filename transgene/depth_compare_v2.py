import logging
import polars as pl
import typer
from pathlib import Path
from typing import Annotated, Optional, List
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from dataclasses import dataclass, asdict

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


@dataclass
class DepthStats:
    sample_id: str
    transgene_depth: float
    background_depth: float
    ratio: float


def read_depth_median(depth_file: Path) -> Optional[float]:
    """
    Reads a depth file and returns the median depth, excluding zero-depth positions.

    Args:
        depth_file: Path to the depth file (expected to be a tab-separated file with a 'Raw Depth' column).

    Returns:
        The median depth as a float, or None if the file cannot be read or processed.
    """
    if not depth_file.exists():
        logger.warning(f"File not found: {depth_file}")
        return None

    try:
        # polars handles gzip compression automatically
        df = pl.read_csv(
            depth_file,
            separator="\t",
            columns=["Raw Depth"] # Only read the necessary column
        )

        depth_series = df["Raw Depth"]

        # Filter out zero depths to consider only covered regions
        valid_depths = depth_series.filter(depth_series > 0)

        if valid_depths.is_empty():
            return 0.0

        return float(valid_depths.median())

    except Exception as e:
        logger.error(f"Error reading {depth_file}: {e}")
        return None


def process_single_sample(
    transgene_dir: Path, background_bamdst_dir: Path
) -> Optional[DepthStats]:
    """
    Processes a single sample to calculate depth statistics.

    Args:
        transgene_dir: Directory containing the transgene sample data.
        background_bamdst_dir: Directory containing the background samples.

    Returns:
        A DepthStats object containing the results, or None if processing fails.
    """
    if not transgene_dir.is_dir():
        return None

    sample_id = transgene_dir.name
    trans_depth_file = transgene_dir / "depth.tsv.gz"
    
    # Assuming the background directory structure matches the sample ID
    bg_depth_file = background_bamdst_dir / sample_id / "depth.tsv.gz"

    transgene_depth = read_depth_median(trans_depth_file)
    
    if transgene_depth is None:
        return None

    if not bg_depth_file.exists():
        logger.debug(f"Background file missing for sample: {sample_id}")
        return None
        
    background_depth = read_depth_median(bg_depth_file)

    if background_depth is None:
        return None

    ratio = (
        transgene_depth / background_depth
        if background_depth > 0
        else 0.0
    )

    return DepthStats(
        sample_id=sample_id,
        transgene_depth=transgene_depth,
        background_depth=background_depth,
        ratio=ratio,
    )


def main(
    trans_gene_bamdst_dir: Annotated[
        Path, typer.Argument(help="Directory containing transgene bamdst results")
    ],
    background_bamdst_dir: Annotated[
        Path, typer.Argument(help="Directory containing background bamdst results")
    ],
    output_file: Annotated[Path, typer.Argument(help="Path to the output Excel file")],
    threads: Annotated[int, typer.Option(help="Number of threads for parallel processing")] = 4,
) -> None:
    """
    Compares the depth of transgene and background BAM files.
    Calculates median depth (excluding 0 coverage) and ratios.
    Outputs the results to an Excel file.
    """

    sample_dirs = [d for d in trans_gene_bamdst_dir.iterdir() if d.is_dir()]
    results: List[DepthStats] = []

    logger.info(
        f"Start processing {len(sample_dirs)} samples with {threads} threads..."
    )

    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(process_single_sample, d, background_bamdst_dir): d.name
            for d in sample_dirs
        }

        for future in tqdm(
            as_completed(futures), total=len(sample_dirs), unit="sample"
        ):
            res = future.result()
            if res:
                results.append(res)

    if not results:
        logger.warning("No results generated.")
        return

    # Sort by sample_id
    results.sort(key=lambda x: x.sample_id)

    # Convert dataclasses to DataFrame
    df_result = pl.DataFrame([asdict(r) for r in results])

    # Ensure output file has .xlsx extension
    if output_file.suffix != ".xlsx":
        output_file = output_file.with_suffix(".xlsx")
        logger.info(f"Changed output extension to .xlsx: {output_file}")

    try:
        df_result.write_excel(output_file, worksheet="Depth Comparison")
        logger.info(f"Results saved to {output_file}")
    except Exception as e:
        logger.error(f"Failed to write Excel file: {e}")
        # Fallback to CSV if Excel fails (e.g., missing engine)
        csv_file = output_file.with_suffix(".tsv")
        df_result.write_csv(csv_file, separator="\t")
        logger.warning(f"Fallback: Results saved to {csv_file}")


if __name__ == "__main__":
    typer.run(main)
