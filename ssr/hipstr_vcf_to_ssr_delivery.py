#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "typer>=0.12",
#   "loguru",
#   "rich",
#   "pandas",
#   "pysam",
#   "primer3-py",
#   "openpyxl",
# ]
# ///
"""Convert HipSTR VCF into polymorphic SSR delivery tables.

This script follows the previously reviewed delivery logic:
- MISA output is not required as input.
- HipSTR VCF is the primary source for locus coordinates, motif metadata,
  diploid repeat genotypes, and VCF-level filtering metrics.
- Primer design and BLAST specificity are optional and only run when the
  required reference resources are supplied.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import pysam
import typer
from loguru import logger
from rich.console import Console
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)

__version__ = "0.1.0"

console = Console()

logger.remove()
logger.add(sys.stderr, level="INFO", format="<level>{level: <8}</level> {message}")

ALLELE_SPLIT_RE = re.compile(r"[|/]")
BLAST_OUTFMT = (
    "6 qseqid sseqid pident length mismatch gapopen qstart qend "
    "sstart send evalue bitscore"
)
BLAST_WORD_SIZE = 11
BLAST_PERC_IDENTITY = 95
BLAST_MAX_TARGET_SEQS = 10
FLANK_SIZE = 400
MIN_TEMPLATE_LENGTH = 100
BASE_COLUMNS = [
    "SSR_ID",
    "seq_id",
    "start",
    "end",
    "motif",
    "unit_size",
    "ref_repeat",
]
STAT_COLUMNS = [
    "maf_like",
    "missing_count",
    "allele_num",
    "repeat_min",
    "repeat_max",
    "repeat_diff",
    "bp_diff",
    "sample_count",
    "het_rate",
]
PRIMER_COLUMNS = [
    "Left_Primer",
    "Right_Primer",
    "Left_Tm",
    "Right_Tm",
    "Product_Size",
    "Specificity",
]


@dataclass
class FilterCriteria:
    alt_alleles_min: int
    alt_alleles_max: int
    pic_min: float
    call_rate_min: float
    skip_initial_filters: bool = False


@dataclass
class ParsedRecord:
    row: dict[str, Any]
    observed_allele_count: int
    pic: float | None
    call_rate: float
    period: int


@dataclass
class RunSummary:
    sample_names: list[str]
    total_input_loci: int = 0
    passed_initial_filters: int = 0
    failed_alt_filter: int = 0
    failed_pic_filter: int = 0
    failed_call_rate_filter: int = 0
    monomorphic_removed: int = 0
    negative_repeat_removed: int = 0
    missing_primer_removed: int = 0
    final_rows: int = 0
    primer_rows_with_sequences: int = 0
    period_counts: Counter[int] = field(default_factory=Counter)
    specificity_counts: Counter[str] = field(default_factory=Counter)


def version_callback(value: bool) -> None:
    if value:
        console.print(f"[bold]hipstr-vcf-to-ssr-delivery[/bold] v{__version__}")
        raise typer.Exit()


def configure_logging(verbose: bool) -> None:
    logger.remove()
    logger.add(
        sys.stderr,
        level="DEBUG" if verbose else "INFO",
        format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{message}</cyan>",
    )


def normalize_fraction(value: float, name: str) -> float:
    if value < 0:
        raise typer.BadParameter(f"{name} must be >= 0")
    if value <= 1:
        return value
    if value <= 100:
        return value / 100
    raise typer.BadParameter(f"{name} must be between 0-1 or 0-100")


def progress_context() -> Progress:
    return Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        console=console,
    )


def derive_motif(ref: str, pos: int, start_info: int, period: int) -> str:
    if not ref or period <= 0:
        return ""
    motif_start = max(0, start_info - pos)
    motif_end = motif_start + period
    motif = ref[motif_start:motif_end]
    if motif:
        return motif
    return ref[:period]


def split_gb_values(gb: Any) -> list[str]:
    if gb is None:
        return []
    return [part for part in ALLELE_SPLIT_RE.split(str(gb)) if part != ""]


def complete_gt(gt: Any) -> tuple[int, ...] | None:
    if gt is None:
        return None
    try:
        alleles = tuple(int(allele) for allele in gt if allele is not None)
    except (TypeError, ValueError):
        return None
    if len(alleles) != len(gt):
        return None
    return alleles


def to_repeat_values(sample_call: Any, ref_repeat: int, period: int) -> list[int] | None:
    gt = sample_call.get("GT")
    gt_complete = complete_gt(gt)
    if gt_complete is None:
        return None

    gb_values = split_gb_values(sample_call.get("GB"))
    if len(gb_values) < len(gt_complete):
        return None
    if any(bp_diff in {".", ""} for bp_diff in gb_values[: len(gt_complete)]):
        return None

    repeats: list[int] = []
    for bp_diff in gb_values[: len(gt_complete)]:
        try:
            bp_diff_int = int(bp_diff)
        except (TypeError, ValueError):
            return None
        repeats.append(ref_repeat + round(bp_diff_int / period))
    return repeats


def format_sample_genotype(repeats: list[int]) -> str:
    if not repeats:
        return ""
    if len(set(repeats)) == 1:
        return str(repeats[0])
    return "/".join(str(value) for value in sorted(repeats))


def calculate_pic(allele_indices: list[int]) -> float | None:
    if not allele_indices:
        return None
    counts = Counter(allele_indices)
    total = sum(counts.values())
    if total == 0:
        return None
    freqs = [count / total for count in counts.values()]
    sum_pi_sq = sum(freq * freq for freq in freqs)
    second_term = 0.0
    for index, freq_i in enumerate(freqs):
        for freq_j in freqs[index + 1 :]:
            second_term += 2 * (freq_i * freq_i) * (freq_j * freq_j)
    pic = 1 - sum_pi_sq - second_term
    return round(pic, 4)


def parse_record(record: Any, sample_names: list[str]) -> ParsedRecord:
    chrom = record.contig
    pos = int(record.pos)
    ref = record.ref or ""
    start_info = int(record.info.get("START", pos))
    end_info = int(record.info.get("END", record.stop))
    period = int(record.info.get("PERIOD", 0))
    if period <= 0:
        raise ValueError(f"Invalid PERIOD for {chrom}:{pos}")

    ssr_id = f"SSR_{chrom}_{pos}"
    locus_end = pos + len(ref) - 1
    motif = derive_motif(ref=ref, pos=pos, start_info=start_info, period=period)
    ref_repeat = max(0, (end_info - start_info + 1) // period)

    all_repeat_alleles: list[int] = []
    sequence_allele_indices: list[int] = []
    sample_values: list[str] = []
    heterozygous_samples = 0

    for sample_name in sample_names:
        sample_call = record.samples[sample_name]
        gt_complete = complete_gt(sample_call.get("GT"))
        if gt_complete is not None:
            sequence_allele_indices.extend(gt_complete)

        repeats = to_repeat_values(sample_call, ref_repeat=ref_repeat, period=period)
        if not repeats:
            sample_values.append("")
            continue

        all_repeat_alleles.extend(repeats)
        if len(set(repeats)) > 1:
            heterozygous_samples += 1
        sample_values.append(format_sample_genotype(repeats))

    sample_count = sum(1 for value in sample_values if value != "")
    missing_count = len(sample_names) - sample_count
    observed_allele_count = len(set(sequence_allele_indices))
    call_rate = sample_count / len(sample_names) if sample_names else 0.0
    pic = calculate_pic(sequence_allele_indices)

    if all_repeat_alleles:
        repeat_counter = Counter(all_repeat_alleles)
        unique_repeats = sorted(repeat_counter)
        repeat_min = unique_repeats[0]
        repeat_max = unique_repeats[-1]
        repeat_diff = repeat_max - repeat_min
        bp_diff = repeat_diff * period
        allele_num = len(unique_repeats)
        maf_like = round(1 - (max(repeat_counter.values()) / len(all_repeat_alleles)), 3)
    else:
        repeat_min = None
        repeat_max = None
        repeat_diff = None
        bp_diff = None
        allele_num = 0
        maf_like = None

    het_rate = round(heterozygous_samples / sample_count, 3) if sample_count > 0 else 0.0

    row = {
        "SSR_ID": ssr_id,
        "seq_id": chrom,
        "start": pos,
        "end": locus_end,
        "motif": motif,
        "unit_size": period,
        "ref_repeat": ref_repeat,
        "maf_like": maf_like,
        "missing_count": missing_count,
        "allele_num": allele_num,
        "repeat_min": repeat_min,
        "repeat_max": repeat_max,
        "repeat_diff": repeat_diff,
        "bp_diff": bp_diff,
        "sample_count": sample_count,
        "het_rate": het_rate,
    }
    row.update(dict(zip(sample_names, sample_values, strict=False)))
    return ParsedRecord(
        row=row,
        observed_allele_count=observed_allele_count,
        pic=pic,
        call_rate=call_rate,
        period=period,
    )


def passes_initial_filters(parsed: ParsedRecord, criteria: FilterCriteria, summary: RunSummary) -> bool:
    if criteria.skip_initial_filters:
        summary.passed_initial_filters += 1
        return True
    fail_alt = not (
        criteria.alt_alleles_min <= parsed.observed_allele_count <= criteria.alt_alleles_max
    )
    fail_pic = parsed.pic is None or parsed.pic < criteria.pic_min
    fail_call_rate = parsed.call_rate < criteria.call_rate_min

    if fail_alt:
        summary.failed_alt_filter += 1
    if fail_pic:
        summary.failed_pic_filter += 1
    if fail_call_rate:
        summary.failed_call_rate_filter += 1

    if fail_alt or fail_pic or fail_call_rate:
        return False

    summary.passed_initial_filters += 1
    return True


def parse_vcf_rows(vcf_path: Path, criteria: FilterCriteria) -> tuple[list[str], list[dict[str, Any]], RunSummary]:
    with pysam.VariantFile(str(vcf_path)) as vcf:
        sample_names = list(vcf.header.samples)
        summary = RunSummary(sample_names=sample_names)
        rows: list[dict[str, Any]] = []

        with progress_context() as progress:
            task = progress.add_task("Parsing HipSTR VCF", total=None)
            for record in vcf:
                summary.total_input_loci += 1
                parsed = parse_record(record, sample_names=sample_names)
                if not passes_initial_filters(parsed, criteria=criteria, summary=summary):
                    progress.advance(task)
                    continue
                if parsed.row["allele_num"] <= 1:
                    summary.monomorphic_removed += 1
                    progress.advance(task)
                    continue
                rows.append(parsed.row)
                progress.advance(task)

    return sample_names, rows, summary


def get_final_columns(sample_names: list[str]) -> list[str]:
    return BASE_COLUMNS + sample_names + STAT_COLUMNS + PRIMER_COLUMNS


def create_dataframe(rows: list[dict[str, Any]], sample_names: list[str]) -> pd.DataFrame:
    final_columns = get_final_columns(sample_names)
    if not rows:
        return pd.DataFrame(columns=final_columns)
    df = pd.DataFrame(rows)
    for column in PRIMER_COLUMNS:
        if column not in df.columns:
            df[column] = ""
    return df


def import_primer3_module() -> Any:
    try:
        import primer3  # type: ignore
    except ImportError as exc:  # pragma: no cover - runtime environment dependent
        raise typer.BadParameter(
            "primer3-py is required when --genome is supplied for primer design."
        ) from exc
    return primer3


def design_primer_for_row(row: pd.Series, fasta: pysam.FastaFile, primer3: Any) -> dict[str, Any]:
    ssr_id = str(row["SSR_ID"])
    chrom = str(row["seq_id"])
    start_1based = int(row["start"])
    end_1based = int(row["end"])
    if start_1based < 1 or end_1based < start_1based:
        return {"SSR_ID": ssr_id, **{column: "" for column in PRIMER_COLUMNS}}

    try:
        chrom_length = fasta.get_reference_length(chrom)
    except KeyError:
        logger.debug(f"Reference contig not found for primer design: {chrom}")
        return {"SSR_ID": ssr_id, **{column: "" for column in PRIMER_COLUMNS}}

    ssr_start0 = start_1based - 1
    ssr_end0 = end_1based
    upstream_start0 = max(0, ssr_start0 - FLANK_SIZE)
    downstream_end0 = min(chrom_length, end_1based + FLANK_SIZE)

    upstream_seq = fasta.fetch(chrom, upstream_start0, ssr_start0)
    ssr_seq = fasta.fetch(chrom, ssr_start0, ssr_end0)
    downstream_seq = fasta.fetch(chrom, end_1based, downstream_end0)

    template_seq = upstream_seq + ssr_seq + downstream_seq
    if len(template_seq) < MIN_TEMPLATE_LENGTH or not ssr_seq:
        return {"SSR_ID": ssr_id, **{column: "" for column in PRIMER_COLUMNS}}

    seq_args = {
        "SEQUENCE_ID": ssr_id,
        "SEQUENCE_TEMPLATE": template_seq,
        "SEQUENCE_TARGET": [len(upstream_seq), len(ssr_seq)],
    }
    global_args = {
        "PRIMER_PRODUCT_SIZE_RANGE": [[100, 300]],
        "PRIMER_OPT_TM": 60.0,
        "PRIMER_MIN_TM": 55.0,
        "PRIMER_MAX_TM": 65.0,
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_MIN_SIZE": 18,
        "PRIMER_MAX_SIZE": 25,
        "PRIMER_NUM_RETURN": 1,
        "PRIMER_MAX_NS_ACCEPTED": 0,
    }

    try:
        if hasattr(primer3.bindings, "designPrimers"):
            result = primer3.bindings.designPrimers(seq_args, global_args)
        else:
            result = primer3.bindings.design_primers(seq_args, global_args)
    except Exception as exc:  # pragma: no cover - primer3 runtime failures are data dependent
        logger.debug(f"Primer3 failed for {ssr_id}: {exc}")
        return {"SSR_ID": ssr_id, **{column: "" for column in PRIMER_COLUMNS}}

    if "PRIMER_LEFT_0_SEQUENCE" not in result or "PRIMER_RIGHT_0_SEQUENCE" not in result:
        return {"SSR_ID": ssr_id, **{column: "" for column in PRIMER_COLUMNS}}

    return {
        "SSR_ID": ssr_id,
        "Left_Primer": str(result.get("PRIMER_LEFT_0_SEQUENCE", "")),
        "Right_Primer": str(result.get("PRIMER_RIGHT_0_SEQUENCE", "")),
        "Left_Tm": round(float(result.get("PRIMER_LEFT_0_TM", 0.0)), 3),
        "Right_Tm": round(float(result.get("PRIMER_RIGHT_0_TM", 0.0)), 3),
        "Product_Size": int(result.get("PRIMER_PAIR_0_PRODUCT_SIZE", 0) or 0),
        "Specificity": "",
    }


def annotate_primers(df: pd.DataFrame, genome: Path) -> pd.DataFrame:
    if df.empty:
        return df

    primer3 = import_primer3_module()
    logger.info(f"Designing primers with genome: {genome}")
    primer_rows: list[dict[str, Any]] = []

    with pysam.FastaFile(str(genome)) as fasta, progress_context() as progress:
        task = progress.add_task("Designing primers", total=len(df))
        for _, row in df.iterrows():
            primer_rows.append(design_primer_for_row(row, fasta=fasta, primer3=primer3))
            progress.advance(task)

    primer_df = pd.DataFrame(primer_rows)
    merged = df.drop(columns=PRIMER_COLUMNS, errors="ignore").merge(primer_df, on="SSR_ID", how="left")
    for column in PRIMER_COLUMNS:
        if column not in merged.columns:
            merged[column] = ""
    merged["Specificity"] = merged["Specificity"].fillna("")
    return merged


def ensure_blast_database(genome: Path, workdir: Path, makeblastdb_bin: str) -> Path:
    if shutil.which(makeblastdb_bin) is None and not Path(makeblastdb_bin).exists():
        raise FileNotFoundError(f"makeblastdb not found: {makeblastdb_bin}")

    blast_dir = workdir / "blastdb"
    blast_dir.mkdir(parents=True, exist_ok=True)
    db_prefix = blast_dir / genome.stem

    existing_extensions = [".nhr", ".nin", ".nsq"]
    if all((db_prefix.with_suffix(ext)).exists() for ext in existing_extensions):
        logger.info(f"Reusing existing BLAST database: {db_prefix}")
        return db_prefix

    logger.info(f"Building BLAST database from genome: {genome}")
    result = subprocess.run(
        [makeblastdb_bin, "-in", str(genome), "-dbtype", "nucl", "-out", str(db_prefix)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        error_text = (result.stderr or result.stdout or "").strip().replace("\n", " | ")
        raise RuntimeError(f"makeblastdb failed with exit code {result.returncode}: {error_text}")
    return db_prefix


def chunked(rows: list[tuple[str, str, str]], size: int) -> list[list[tuple[str, str, str]]]:
    return [rows[index : index + size] for index in range(0, len(rows), size)]


def run_blast_batch(
    rows: list[tuple[str, str, str]],
    blastn_bin: str,
    blast_db: Path,
) -> dict[str, str]:
    if shutil.which(blastn_bin) is None and not Path(blastn_bin).exists():
        raise FileNotFoundError(f"blastn not found: {blastn_bin}")

    fasta_path: Path | None = None
    try:
        query_lengths: dict[str, int] = {}
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", suffix=".fasta", delete=False) as handle:
            fasta_path = Path(handle.name)
            for ssr_id, left_primer, right_primer in rows:
                left_id = f"{ssr_id}_L"
                right_id = f"{ssr_id}_R"
                query_lengths[left_id] = len(left_primer)
                query_lengths[right_id] = len(right_primer)
                handle.write(f">{left_id}\n{left_primer}\n")
                handle.write(f">{right_id}\n{right_primer}\n")

        result = subprocess.run(
            [
                blastn_bin,
                "-query",
                str(fasta_path),
                "-db",
                str(blast_db),
                "-task",
                "blastn-short",
                "-evalue",
                "1",
                "-perc_identity",
                str(BLAST_PERC_IDENTITY),
                "-word_size",
                str(BLAST_WORD_SIZE),
                "-max_target_seqs",
                str(BLAST_MAX_TARGET_SEQS),
                "-outfmt",
                BLAST_OUTFMT,
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            error_text = (result.stderr or result.stdout or "").strip().replace("\n", " | ")
            raise RuntimeError(f"blastn failed with exit code {result.returncode}: {error_text}")

        perfect_hits: dict[str, int] = {query_id: 0 for query_id in query_lengths}
        for line in result.stdout.splitlines():
            if not line:
                continue
            columns = line.split("\t")
            if len(columns) < 4:
                continue
            query_id = columns[0]
            query_length = query_lengths.get(query_id)
            if query_length is None:
                continue
            try:
                pident = float(columns[2])
                alignment_length = int(columns[3])
            except ValueError:
                continue
            min_alignment_length = max(query_length - 1, round(query_length * 0.95))
            if pident >= BLAST_PERC_IDENTITY and alignment_length >= min_alignment_length:
                perfect_hits[query_id] += 1

        batch_specificity: dict[str, str] = {}
        for ssr_id, _, _ in rows:
            left_hits = perfect_hits.get(f"{ssr_id}_L", 0)
            right_hits = perfect_hits.get(f"{ssr_id}_R", 0)
            if left_hits == 1 and right_hits == 1:
                batch_specificity[ssr_id] = "Specific"
            else:
                batch_specificity[ssr_id] = f"Non-specific (L:{left_hits}, R:{right_hits})"
        return batch_specificity
    finally:
        if fasta_path is not None:
            fasta_path.unlink(missing_ok=True)


def annotate_specificity(
    df: pd.DataFrame,
    blast_db: Path,
    blastn_bin: str,
    batch_size: int,
) -> pd.DataFrame:
    if df.empty:
        return df

    primer_rows: list[tuple[str, str, str]] = []
    for _, row in df.iterrows():
        left_primer = str(row.get("Left_Primer", "") or "")
        right_primer = str(row.get("Right_Primer", "") or "")
        if left_primer and right_primer:
            primer_rows.append((str(row["SSR_ID"]), left_primer, right_primer))

    if not primer_rows:
        logger.warning("No designed primers available; skipping specificity annotation.")
        return df

    logger.info(f"Annotating primer specificity with BLAST DB: {blast_db}")
    specificity_map: dict[str, str] = {}
    batches = chunked(primer_rows, size=batch_size)

    with progress_context() as progress:
        task = progress.add_task("Checking specificity", total=len(batches))
        for batch in batches:
            specificity_map.update(run_blast_batch(batch, blastn_bin=blastn_bin, blast_db=blast_db))
            progress.advance(task)

    annotated = df.copy()
    annotated["Specificity"] = annotated["SSR_ID"].map(specificity_map).fillna("")
    return annotated


def clean_final_dataframe(
    df: pd.DataFrame,
    summary: RunSummary,
    clean_final: bool,
    primers_attempted: bool,
) -> pd.DataFrame:
    if not clean_final or df.empty:
        return df

    cleaned = df.copy()
    initial_count = len(cleaned)

    if "repeat_min" in cleaned.columns:
        repeat_min_numeric = pd.to_numeric(cleaned["repeat_min"], errors="coerce")
        cleaned = cleaned[repeat_min_numeric.isna() | (repeat_min_numeric >= 0)].copy()
    summary.negative_repeat_removed = initial_count - len(cleaned)

    if primers_attempted:
        before_primers = len(cleaned)
        left_primer = cleaned["Left_Primer"].fillna("").astype(str)
        right_primer = cleaned["Right_Primer"].fillna("").astype(str)
        cleaned = cleaned[(left_primer != "") & (right_primer != "")].copy()
        summary.missing_primer_removed = before_primers - len(cleaned)

    return cleaned.reset_index(drop=True)


def format_output_dataframe(df: pd.DataFrame, sample_names: list[str]) -> pd.DataFrame:
    output = df.copy()
    final_columns = get_final_columns(sample_names)
    for column in final_columns:
        if column not in output.columns:
            output[column] = ""

    numeric_round3 = ["maf_like", "het_rate", "Left_Tm", "Right_Tm"]
    integer_columns = [
        "start",
        "end",
        "unit_size",
        "ref_repeat",
        "missing_count",
        "allele_num",
        "repeat_min",
        "repeat_max",
        "repeat_diff",
        "bp_diff",
        "sample_count",
        "Product_Size",
    ]

    for column in numeric_round3:
        output[column] = pd.to_numeric(output[column], errors="coerce").round(3)

    for column in integer_columns:
        output[column] = pd.to_numeric(output[column], errors="coerce").astype("Int64")

    output = output[final_columns].copy()
    output = output.astype(object).where(pd.notna(output), "")
    return output


def write_outputs(df: pd.DataFrame, out_tsv: Path, out_xlsx: Path) -> None:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False)
    df.to_excel(out_xlsx, index=False)


def update_summary_from_output(df: pd.DataFrame, summary: RunSummary) -> None:
    summary.final_rows = len(df)
    if "unit_size" in df.columns:
        summary.period_counts = Counter(
            int(value)
            for value in pd.to_numeric(df["unit_size"], errors="coerce").dropna().astype(int).tolist()
        )
    if "Left_Primer" in df.columns:
        summary.primer_rows_with_sequences = int((df["Left_Primer"].fillna("").astype(str) != "").sum())
    if "Specificity" in df.columns:
        labels = [str(value).strip() for value in df["Specificity"].fillna("").tolist() if str(value).strip()]
        summary.specificity_counts = Counter(labels)


def render_report(
    vcf_path: Path,
    report_path: Path,
    criteria: FilterCriteria,
    summary: RunSummary,
    genome: Path | None,
    out_tsv: Path,
    out_xlsx: Path,
) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)

    if criteria.skip_initial_filters:
        filtering_lines = [
            "Filtering Criteria:",
            "  - Initial VCF-level filtering: skipped",
        ]
    else:
        filtering_lines = [
            "Filtering Criteria:",
            (
                "  - Alleles per locus "
                f"(observed sequence alleles from GT): {criteria.alt_alleles_min}-{criteria.alt_alleles_max}"
            ),
            f"  - PIC >= {criteria.pic_min:.2f}",
            f"  - Call rate >= {criteria.call_rate_min * 100:.2f}%",
        ]

    period_lines = []
    if summary.period_counts:
        period_lines.append("By Repeat Unit Length:")
        for period, count in sorted(summary.period_counts.items()):
            period_lines.append(f"  - Period {period}: {count} loci")

    specificity_lines = []
    if summary.specificity_counts:
        specificity_lines.append("Primer Specificity:")
        for label, count in summary.specificity_counts.items():
            specificity_lines.append(f"  - {label}: {count}")

    content = "\n".join(
        [
            "SSR Filtering Report",
            "=" * 50,
            "",
            f"Input VCF: {vcf_path}",
            f"Total input loci: {summary.total_input_loci}",
            f"Samples: {', '.join(summary.sample_names)}",
            "",
            *filtering_lines,
            "",
            "Filtering Results:",
            f"  - Passed initial filters: {summary.passed_initial_filters} loci",
            f"  - Failed ALT allele filter: {summary.failed_alt_filter}",
            f"  - Failed PIC filter: {summary.failed_pic_filter}",
            f"  - Failed call rate filter: {summary.failed_call_rate_filter}",
            "  - Note: failure counts above are not mutually exclusive.",
            "",
            "Post-processing Notes",
            "=" * 50,
            "",
            "Monomorphic site removal:",
            f"  - {summary.monomorphic_removed} loci with allele_num <= 1 were removed.",
            "",
            "Final cleaning:",
            f"  - repeat_min < 0 removed: {summary.negative_repeat_removed}",
            f"  - missing primers removed: {summary.missing_primer_removed}",
            "",
            "Final output:",
            f"  - Rows: {summary.final_rows}",
            f"  - TSV: {out_tsv}",
            f"  - XLSX: {out_xlsx}",
            f"  - Genome used for primers: {genome if genome else 'No'}",
            f"  - Rows with designed primers: {summary.primer_rows_with_sequences}",
            "",
            *period_lines,
            "",
            *specificity_lines,
            "",
        ]
    )
    report_path.write_text(content, encoding="utf-8")


def run_pipeline(
    vcf: Path,
    out_xlsx: Path,
    genome: Path | None,
    blast_db: Path | None,
    workdir: Path | None,
    out_tsv: Path | None,
    report_txt: Path | None,
    alleles_min: int,
    alleles_max: int,
    pic_min: float,
    call_rate_min: float,
    skip_initial_filters: bool,
    clean_final: bool,
    check_specificity: bool,
    batch_size: int,
    blastn: str,
    makeblastdb: str,
) -> None:
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    resolved_out_tsv = out_tsv or out_xlsx.with_suffix(".tsv")
    resolved_report_txt = report_txt or out_xlsx.with_name("filtering_report.txt")
    normalized_call_rate = normalize_fraction(call_rate_min, "call_rate_min")

    criteria = FilterCriteria(
        alt_alleles_min=alleles_min,
        alt_alleles_max=alleles_max,
        pic_min=pic_min,
        call_rate_min=normalized_call_rate,
        skip_initial_filters=skip_initial_filters,
    )

    logger.info(f"Input VCF: {vcf}")
    sample_names, rows, summary = parse_vcf_rows(vcf_path=vcf, criteria=criteria)
    logger.info(f"Rows retained after VCF parsing: {len(rows)}")

    df = create_dataframe(rows=rows, sample_names=sample_names)

    primers_attempted = False
    if genome is not None:
        primers_attempted = True
        df = annotate_primers(df=df, genome=genome)
    else:
        logger.info("No genome supplied; primer columns will remain empty.")

    if check_specificity and not df.empty:
        specificity_source = blast_db
        if specificity_source is None and genome is not None:
            resolved_workdir = workdir or out_xlsx.parent / ".ssr_delivery_work"
            specificity_source = ensure_blast_database(
                genome=genome,
                workdir=resolved_workdir,
                makeblastdb_bin=makeblastdb,
            )

        if specificity_source is not None:
            try:
                df = annotate_specificity(
                    df=df,
                    blast_db=specificity_source,
                    blastn_bin=blastn,
                    batch_size=batch_size,
                )
            except Exception as exc:
                logger.warning(f"Specificity annotation skipped: {exc}")
        else:
            logger.info("No BLAST DB or genome supplied; specificity annotation skipped.")

    df = clean_final_dataframe(
        df=df,
        summary=summary,
        clean_final=clean_final,
        primers_attempted=primers_attempted,
    )
    df = format_output_dataframe(df=df, sample_names=sample_names)
    update_summary_from_output(df=df, summary=summary)

    write_outputs(df=df, out_tsv=resolved_out_tsv, out_xlsx=out_xlsx)
    render_report(
        vcf_path=vcf,
        report_path=resolved_report_txt,
        criteria=criteria,
        summary=summary,
        genome=genome,
        out_tsv=resolved_out_tsv,
        out_xlsx=out_xlsx,
    )

    console.print("\n[green]✓ Delivery table generated[/green]")
    console.print(f"  Rows: {summary.final_rows}")
    console.print(f"  Samples: {len(sample_names)}")
    console.print(f"  TSV: {resolved_out_tsv}")
    console.print(f"  XLSX: {out_xlsx}")
    console.print(f"  Report: {resolved_report_txt}")


def main(
    vcf: Annotated[
        Path,
        typer.Argument(help="Input HipSTR VCF(.vcf or .vcf.gz)", exists=True, dir_okay=False),
    ],
    out_xlsx: Annotated[
        Path,
        typer.Argument(help="Output xlsx path, e.g. polymorphic_SSRs_with_primers.maf_like.allele_num.xlsx"),
    ],
    genome: Annotated[
        Path | None,
        typer.Option(
            "--genome",
            help="Reference genome FASTA used for primer design. If omitted, primer columns stay empty.",
        ),
    ] = None,
    blast_db: Annotated[
        Path | None,
        typer.Option("--blast-db", help="Existing BLAST database prefix for primer specificity annotation."),
    ] = None,
    workdir: Annotated[
        Path | None,
        typer.Option("--workdir", help="Optional work/cache directory for auto-generated BLAST databases."),
    ] = None,
    out_tsv: Annotated[
        Path | None,
        typer.Option("--out-tsv", help="Optional TSV output path. Defaults to the xlsx sibling with .tsv suffix."),
    ] = None,
    report_txt: Annotated[
        Path | None,
        typer.Option("--report-txt", help="Optional filtering report path. Defaults to filtering_report.txt beside the xlsx."),
    ] = None,
    alleles_min: Annotated[
        int,
        typer.Option(
            "--alleles-min",
            help="Minimum observed sequence allele count from GT for initial filtering.",
        ),
    ] = 4,
    alleles_max: Annotated[
        int,
        typer.Option(
            "--alleles-max",
            help="Maximum observed sequence allele count from GT for initial filtering.",
        ),
    ] = 7,
    pic_min: Annotated[
        float,
        typer.Option("--pic-min", help="Minimum PIC threshold for initial filtering."),
    ] = 0.25,
    call_rate_min: Annotated[
        float,
        typer.Option("--call-rate-min", help="Minimum call rate threshold; accepts 0.8 or 80."),
    ] = 0.8,
    skip_initial_filters: Annotated[
        bool,
        typer.Option(
            "--skip-initial-filters",
            help="Skip observed-allele-count / PIC / call-rate filtering.",
        ),
    ] = False,
    clean_final: Annotated[
        bool,
        typer.Option(
            "--clean-final/--no-clean-final",
            help="Remove repeat_min < 0; when primers are designed, also drop rows without primer pairs.",
        ),
    ] = True,
    check_specificity: Annotated[
        bool,
        typer.Option(
            "--check-specificity/--no-check-specificity",
            help="Attempt BLAST specificity annotation when a BLAST DB or genome is available.",
        ),
    ] = True,
    batch_size: Annotated[
        int,
        typer.Option("--batch-size", help="Number of primer pairs per BLAST batch."),
    ] = 100,
    blastn: Annotated[
        str,
        typer.Option("--blastn", help="blastn executable path."),
    ] = "blastn",
    makeblastdb: Annotated[
        str,
        typer.Option("--makeblastdb", help="makeblastdb executable path."),
    ] = "makeblastdb",
    verbose: Annotated[
        bool,
        typer.Option("--verbose", "-v", help="Enable debug logging."),
    ] = False,
    version: Annotated[
        bool | None,
        typer.Option("--version", callback=version_callback, is_eager=True, help="Show version and exit."),
    ] = None,
) -> None:
    """Convert HipSTR VCF into delivery TSV/XLSX tables."""
    configure_logging(verbose=verbose)

    if genome is not None and not genome.exists():
        raise typer.BadParameter(f"Genome FASTA not found: {genome}")
    if batch_size <= 0:
        raise typer.BadParameter("--batch-size must be a positive integer")
    if alleles_min > alleles_max:
        raise typer.BadParameter("--alleles-min cannot be greater than --alleles-max")

    try:
        run_pipeline(
            vcf=vcf,
            out_xlsx=out_xlsx,
            genome=genome,
            blast_db=blast_db,
            workdir=workdir,
            out_tsv=out_tsv,
            report_txt=report_txt,
            alleles_min=alleles_min,
            alleles_max=alleles_max,
            pic_min=pic_min,
            call_rate_min=call_rate_min,
            skip_initial_filters=skip_initial_filters,
            clean_final=clean_final,
            check_specificity=check_specificity,
            batch_size=batch_size,
            blastn=blastn,
            makeblastdb=makeblastdb,
        )
    except typer.BadParameter:
        raise
    except Exception as exc:
        logger.error(f"Pipeline failed: {exc}")
        raise typer.Exit(code=1) from exc


if __name__ == "__main__":
    typer.run(main)
