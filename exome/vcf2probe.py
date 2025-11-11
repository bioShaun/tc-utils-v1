"""Convert SNP records to probe sequences with variant/reference contexts."""

from dataclasses import dataclass
from pathlib import Path
from contextlib import ExitStack
from typing import Iterator, Optional, TextIO, Tuple

import pandas as pd
import typer
from pyfaidx import Fasta
from tqdm import tqdm

CHUNK_SIZE = 1_000_000
TABLE_NAME = "seq.table.csv"
VARIANT_FASTA_NAME = "seq.fasta"
REFERENCE_FASTA_NAME = "seq.reference.fasta"


VALID_SEQUENCE_MODES = ("variant", "reference", "both")


@dataclass
class SequenceOutputs:
    """Paths to the generated table and FASTA artifacts."""

    table: Path
    variant_fasta: Optional[Path]
    reference_fasta: Optional[Path]


def build_sequence_strings(
    row: pd.Series, reference: Fasta, half_length: int
) -> Tuple[str, str]:
    """
    Return (variant_sequence, reference_sequence) around a SNP.
    """
    pos = int(row["pos"]) - 1
    chrom = str(row["chrom"])
    if chrom not in reference:
        raise KeyError(f"Chromosome {chrom} 不在参考序列中")

    chrom_seq = reference[chrom]
    chrom_len = len(chrom_seq)
    if pos < 0 or pos >= chrom_len:
        raise ValueError(f"位点 {chrom}:{pos + 1} 超出参考序列范围 (长度 {chrom_len})")

    start = max(0, pos - half_length)
    end = min(chrom_len, pos + 1 + half_length)

    left_seq = chrom_seq[start:pos].seq if start < pos else ""
    right_seq = chrom_seq[pos + 1 : end].seq if (pos + 1) < end else ""

    right_seq = right_seq or ""
    variant_seq = f"{left_seq}[{row['ref']}/{row['alt']}]"
    reference_seq = f"{left_seq}{row['ref']}"
    return f"{variant_seq}{right_seq}", f"{reference_seq}{right_seq}"


def normalize_sequence_mode(sequence_mode: str) -> str:
    """
    Validate and normalize the requested sequence mode.
    """
    normalized = sequence_mode.lower()
    if normalized not in VALID_SEQUENCE_MODES:
        raise ValueError(
            f"sequence_mode 只支持: {', '.join(VALID_SEQUENCE_MODES)}; 当前为 {sequence_mode}"
        )
    return normalized


def maybe_with_progress(iterable: Iterator[pd.DataFrame], enable: bool):
    """
    Wrap iterator with tqdm when progress reporting is enabled.
    """
    if enable:
        return tqdm(iterable, desc="Processing chunks", unit="chunk")
    return iterable


def iter_variant_chunks(vcf: Path, is_vcf: bool) -> Iterator[pd.DataFrame]:
    """
    Yield chunks of variant rows from either VCF or TSV input.
    """
    if is_vcf:
        reader = pd.read_table(
            vcf,
            chunksize=CHUNK_SIZE,
            header=None,
            usecols=[0, 1, 3, 4],
            names=["chrom", "pos", "ref", "alt"],
            comment="#",
        )
    else:
        reader = pd.read_csv(
            vcf,
            chunksize=CHUNK_SIZE,
            header=None,
            names=["chrom", "pos", "ref", "alt"],
        )
    yield from reader


def filter_snps(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep only SNP rows.
    """
    mask = (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)
    return df.loc[mask].copy()


def prepare_output_files(out_dir: Path, sequence_mode: str) -> SequenceOutputs:
    """
    Ensure the output directory exists and previous run artifacts are removed.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    variant_path = out_dir / VARIANT_FASTA_NAME if sequence_mode in {"variant", "both"} else None
    reference_path = out_dir / REFERENCE_FASTA_NAME if sequence_mode in {"reference", "both"} else None
    outputs = SequenceOutputs(
        table=out_dir / TABLE_NAME,
        variant_fasta=variant_path,
        reference_fasta=reference_path,
    )
    for path in (outputs.table, outputs.variant_fasta, outputs.reference_fasta):
        if path and path.exists():
            path.unlink()
    return outputs


def write_chunk_outputs(
    df: pd.DataFrame,
    outputs: SequenceOutputs,
    variant_handle: Optional[TextIO],
    reference_handle: Optional[TextIO],
    include_header: bool,
) -> None:
    """
    Append the chunk to the summary table and FASTA files.
    """

    df.to_csv(outputs.table, header=include_header, index=False, mode="a")
    if variant_handle is not None:
        for chrom, pos, sequence in zip(df["chrom"], df["pos"], df["sequence"]):
            seq_id = f"{chrom}_{pos}"
            variant_handle.write(f">{seq_id}\n{sequence}\n")
    if reference_handle is not None:
        for chrom, pos, reference_sequence in zip(
            df["chrom"], df["pos"], df["reference_sequence"]
        ):
            seq_id = f"{chrom}_{pos}"
            reference_handle.write(f">{seq_id}\n{reference_sequence}\n")


def run(
    vcf: Path,
    ref: Path,
    out_dir: Path,
    half_length: int = 200,
    is_vcf: bool = True,
    sequence_mode: str = "both",
    show_progress: bool = True,
) -> None:
    """
    Stream variants and build FASTA/table outputs that contain variant and/or
    reference context sequences depending on `sequence_mode`.
    """
    normalized_mode = normalize_sequence_mode(sequence_mode)
    reference = Fasta(str(ref))
    outputs = prepare_output_files(out_dir, normalized_mode)
    table_written = False

    try:
        with ExitStack() as stack:
            variant_handle = (
                stack.enter_context(outputs.variant_fasta.open("a", encoding="utf-8"))
                if outputs.variant_fasta
                else None
            )
            reference_handle = (
                stack.enter_context(
                    outputs.reference_fasta.open("a", encoding="utf-8")
                )
                if outputs.reference_fasta
                else None
            )
            chunk_iter = maybe_with_progress(iter_variant_chunks(vcf, is_vcf), show_progress)
            for chunk in chunk_iter:
                snp_df = filter_snps(chunk)
                if snp_df.empty:
                    continue
                snp_df[["sequence", "reference_sequence"]] = snp_df.apply(
                    lambda row: build_sequence_strings(row, reference, half_length),
                    axis=1,
                    result_type="expand",
                )
                write_chunk_outputs(
                    snp_df,
                    outputs,
                    variant_handle,
                    reference_handle,
                    include_header=not table_written,
                )
                table_written = True
    finally:
        reference.close()


def main(
    vcf: Path,
    ref: Path,
    out_dir: Path,
    half_length: int = 200,
    is_vcf: bool = True,
    sequence_mode: str = typer.Option(
        "both",
        "--sequence-mode",
        "-m",
        help="选择输出 variant、reference 或 both 序列",
    ),
    show_progress: bool = typer.Option(
        True,
        "--show-progress/--no-progress",
        help="是否使用 tqdm 显示进度条",
    ),
) -> None:
    """
    CLI entry point compatible with typer.
    """
    run(vcf, ref, out_dir, half_length, is_vcf, sequence_mode, show_progress)


if __name__ == "__main__":
    typer.run(main)
