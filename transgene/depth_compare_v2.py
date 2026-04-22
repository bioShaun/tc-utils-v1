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
    transgene_coverage: float
    background_depth: float
    ratio: float
    genotype: str


def classify_genotype(
    transgene_coverage: float,
    ratio: float,
    tolerance: float,
    min_coverage: float = 0.3,
    het_ratio: float = 0.75,
    hom_ratio: float = 1.0,
) -> str:
    """根据 transgene_coverage 和 ratio 判断转基因纯合/杂合/非转基因。

    - coverage < min_coverage: 非转基因
    - coverage >= min_coverage 且 ratio 在 het_ratio ± tolerance: 杂合
    - coverage >= min_coverage 且 ratio >= hom_ratio: 纯合（上不封顶）
    """
    if transgene_coverage < min_coverage:
        return "非转基因"
    if abs(ratio - het_ratio) <= tolerance:
        return "转基因杂合"
    if ratio >= hom_ratio:
        return "转基因纯合"
    return "未确定"


@dataclass
class DepthResult:
    """read_depth_file 的返回结果。"""
    median: float
    coverage: float


def read_depth_file(depth_file: Path) -> Optional[DepthResult]:
    """
    Reads a depth file and returns the median depth (excluding zero-depth positions)
    and coverage ratio (proportion of positions with depth > 0).

    Args:
        depth_file: Path to the depth file (tab-separated with a 'Raw Depth' column).

    Returns:
        A DepthResult with median and coverage, or None if the file cannot be read.
    """
    if not depth_file.exists():
        logger.warning(f"File not found: {depth_file}")
        return None

    try:
        # polars handles gzip compression automatically
        df = pl.read_csv(
            depth_file,
            separator="\t",
            columns=["Raw Depth"],
        )

        depth_series = df["Raw Depth"]
        total = len(depth_series)
        if total == 0:
            return DepthResult(median=0.0, coverage=0.0)

        valid_depths = depth_series.filter(depth_series > 0)
        coverage = len(valid_depths) / total
        median = float(valid_depths.median()) if not valid_depths.is_empty() else 0.0

        return DepthResult(median=median, coverage=coverage)

    except Exception as e:
        logger.error(f"Error reading {depth_file}: {e}")
        return None


def process_single_sample(
    transgene_dir: Path,
    background_bamdst_dir: Path,
    tolerance: float,
    min_coverage: float,
    het_ratio: float,
    hom_ratio: float,
) -> Optional[DepthStats]:
    """
    Processes a single sample to calculate depth statistics.

    Args:
        transgene_dir: Directory containing the transgene sample data.
        background_bamdst_dir: Directory containing the background samples.
        tolerance: Tolerance for genotype classification.

    Returns:
        A DepthStats object containing the results, or None if processing fails.
    """
    if not transgene_dir.is_dir():
        return None

    sample_id = transgene_dir.name
    trans_depth_file = transgene_dir / "depth.tsv.gz"

    # Assuming the background directory structure matches the sample ID
    bg_depth_file = background_bamdst_dir / sample_id / "depth.tsv.gz"

    trans_result = read_depth_file(trans_depth_file)

    if trans_result is None:
        return None

    if not bg_depth_file.exists():
        logger.debug(f"Background file missing for sample: {sample_id}")
        return None

    bg_result = read_depth_file(bg_depth_file)

    if bg_result is None:
        return None

    ratio = (
        trans_result.median / bg_result.median if bg_result.median > 0 else 0.0
    )

    genotype = classify_genotype(
        trans_result.coverage, ratio, tolerance, min_coverage, het_ratio, hom_ratio
    )

    return DepthStats(
        sample_id=sample_id,
        transgene_depth=trans_result.median,
        transgene_coverage=trans_result.coverage,
        background_depth=bg_result.median,
        ratio=ratio,
        genotype=genotype,
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
    tolerance: Annotated[
        float, typer.Option(help="杂合判定的浮动容差范围")
    ] = 0.15,
    min_coverage: Annotated[
        float, typer.Option(help="判定为转基因所需的最低 transgene coverage")
    ] = 0.3,
    het_ratio: Annotated[
        float, typer.Option(help="转基因杂合的 ratio 中心值")
    ] = 0.75,
    hom_ratio: Annotated[
        float, typer.Option(help="转基因纯合的 ratio 下限阈值（>=）")
    ] = 1.0,
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
            executor.submit(
                process_single_sample, d, background_bamdst_dir,
                tolerance, min_coverage, het_ratio, hom_ratio,
            ): d.name
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
