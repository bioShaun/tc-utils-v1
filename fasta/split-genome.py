from __future__ import annotations

import textwrap
from typing import Iterable, Tuple
from pathlib import Path

import pandas as pd
import typer
from loguru import logger
from pyfaidx import Fasta
from tqdm import tqdm

SPLIT_SIZE = 500_000_000
FASTA_LINE_LENGTH = 60
PRIORITY_FEATURES: Tuple[str, ...] = ("gene", "mRNA", "transcript", "CDS", "exon")

GFF_COLUMNS = [
    "seqid",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
]
BED_COLUMNS = ["Chromosome", "start", "end", "id"]


def parse_gff_with_pandas(
    gff_file: Path | str,
    feature_type: str | None = None,
    priority_features: Iterable[str] = PRIORITY_FEATURES,
) -> pd.DataFrame:
    """Load a GFF file and return records for the selected feature type."""
    df = pd.read_csv(
        gff_file,
        sep="\t",
        comment="#",
        names=GFF_COLUMNS,
        na_values=".",
        skip_blank_lines=True,
    )
    if df.empty:
        raise ValueError("GFF file contains no feature records.")

    if feature_type:
        genes_df = df[df["feature"] == feature_type].copy()
        if genes_df.empty:
            raise ValueError(f"Feature type '{feature_type}' not found in GFF file.")
        logger.info("使用指定的feature类型: {} (共 {} 个)", feature_type, len(genes_df))
    else:
        feature_counts = df["feature"].value_counts()
        selected_feature = None
        for feat in priority_features:
            if feat in feature_counts.index and feature_counts[feat] > 0:
                selected_feature = feat
                break
        if selected_feature is None:
            selected_feature = feature_counts.index[0]
        genes_df = df[df["feature"] == selected_feature].copy()
        logger.info("自动选择feature类型: {} (共 {} 个)", selected_feature, len(genes_df))
        logger.debug(
            "feature类型Top10: {}",
            ", ".join(feature_counts.index[:10].tolist()),
        )

    genes_df["start"] = genes_df["start"].astype(int)
    genes_df["end"] = genes_df["end"].astype(int)
    return genes_df


def calculate_gaps(genes_df: pd.DataFrame, show_progress: bool = True) -> pd.DataFrame:
    """Calculate intergenic gaps."""
    if genes_df.empty:
        return pd.DataFrame(
            columns=[
                "Chromosome",
                "gap_size",
                "gap_start",
                "gap_end",
                "upstream_gene_end",
                "downstream_gene_start",
                "upstream_gene_strand",
                "downstream_gene_strand",
                "upstream_gene_info",
                "downstream_gene_info",
            ]
        )

    genes_sorted = genes_df.sort_values(["seqid", "start"]).reset_index(drop=True)
    gaps = []
    for seqid, group in genes_sorted.groupby("seqid"):
        group = group.reset_index(drop=True)
        iterator = range(len(group) - 1)
        if show_progress:
            iterator = tqdm(iterator, desc=f"计算染色体{seqid}的基因间隔")

        for i in iterator:
            gene1 = group.iloc[i]
            gene2 = group.iloc[i + 1]
            gap_size = gene2["start"] - gene1["end"] - 1
            if gap_size > 0:
                gaps.append(
                    {
                        "Chromosome": seqid,
                        "gap_size": gap_size,
                        "gap_start": gene1["end"] + 1,
                        "gap_end": gene2["start"] - 1,
                        "upstream_gene_end": gene1["end"],
                        "downstream_gene_start": gene2["start"],
                        "upstream_gene_strand": gene1["strand"],
                        "downstream_gene_strand": gene2["strand"],
                        "upstream_gene_info": gene1["attributes"],
                        "downstream_gene_info": gene2["attributes"],
                    }
                )

    return pd.DataFrame(gaps)


def generate_split_chr_bed(best_split_site_df: pd.DataFrame) -> pd.DataFrame:
    """Convert best gap locations into BED-like split definitions."""
    if best_split_site_df.empty:
        return pd.DataFrame(columns=BED_COLUMNS)

    records = []
    for row in best_split_site_df.itertuples(index=False):
        split_point = (row.gap_start + row.gap_end) // 2
        records.extend(
            [
                {
                    "Chromosome": row.Chromosome,
                    "start": 0,
                    "end": split_point,
                    "id": f"{row.Chromosome}a",
                },
                {
                    "Chromosome": row.Chromosome,
                    "start": split_point,
                    "end": row.chrom_size,
                    "id": f"{row.Chromosome}b",
                },
            ]
        )
    return pd.DataFrame(records, columns=BED_COLUMNS)


def select_split_candidates(
    chrom_bounds: pd.DataFrame,
    gap_df: pd.DataFrame,
    min_gene_gap: int,
) -> pd.DataFrame:
    """Select the best gap per chromosome that satisfies the minimum gap size."""
    if chrom_bounds.empty:
        return pd.DataFrame(
            columns=[
                "Chromosome",
                "gap_size",
                "gap_start",
                "gap_end",
                "chrom_size",
                "split_start",
                "split_end",
            ]
        )

    merged = gap_df.merge(chrom_bounds, on="Chromosome", how="inner")
    window_mask = (merged["gap_start"] >= merged["split_start"]) & (
        merged["gap_end"] <= merged["split_end"]
    )
    candidates = merged[window_mask].copy()
    if candidates.empty:
        missing = ", ".join(sorted(chrom_bounds["Chromosome"].unique()))
        raise ValueError(
            f"No suitable gaps found within split window for chromosomes: {missing}"
        )

    best_idx = candidates.groupby("Chromosome")["gap_size"].idxmax()
    best = candidates.loc[best_idx].copy()

    missing_chromosomes = set(chrom_bounds["Chromosome"]) - set(best["Chromosome"])
    if missing_chromosomes:
        raise ValueError(
            "无法为以下染色体找到合适的切割位置: "
            f"{', '.join(sorted(missing_chromosomes))}"
        )

    too_small = best[best["gap_size"] < min_gene_gap]
    if not too_small.empty:
        names = ", ".join(sorted(too_small["Chromosome"].unique()))
        raise ValueError(
            "以下染色体的基因间隔小于最小阈值(min_gene_gap): "
            f"{names} (阈值: {min_gene_gap})"
        )
    return best


def generate_split_gff(split_chr_bed: pd.DataFrame, gff: Path, out_gff: Path) -> None:
    """Write a split GFF file by remapping coordinates to the new contigs."""
    if split_chr_bed.empty:
        logger.warning("没有需要切割的染色体，跳过生成拆分后的 GFF 文件。")
        return

    gff_df = pd.read_table(
        gff,
        sep="\t",
        header=None,
        names=GFF_COLUMNS,
        skip_blank_lines=True,
        comment="#",
    )
    rename_split_chr_bed = split_chr_bed.rename(
        columns={"Chromosome": "seqid", "start": "split_start", "end": "split_end", "id": "new_seqid"}
    )
    add_split_chr_df = gff_df.merge(rename_split_chr_bed, on="seqid", how="inner")
    filter1 = add_split_chr_df["start"] >= add_split_chr_df["split_start"]
    filter2 = add_split_chr_df["start"] <= add_split_chr_df["split_end"]
    filter_add_split_chr_df = add_split_chr_df[filter1 & filter2].copy()
    filter_add_split_chr_df["new_start"] = (
        filter_add_split_chr_df["start"] - filter_add_split_chr_df["split_start"]
    )
    filter_add_split_chr_df["new_end"] = (
        filter_add_split_chr_df["end"] - filter_add_split_chr_df["split_start"]
    )
    filter_add_split_chr_df.sort_values(
        ["new_seqid", "new_start", "new_end"], inplace=True
    )
    filter_add_split_chr_df.to_csv(
        out_gff,
        sep="\t",
        index=False,
        header=False,
        columns=[
            "new_seqid",
            "source",
            "feature",
            "new_start",
            "new_end",
            "score",
            "strand",
            "phase",
            "attributes",
        ],
    )


def generate_split_genome(
    fasta_path: Path,
    bed_df: pd.DataFrame,
    out_fasta_path: Path,
    line_length: int = FASTA_LINE_LENGTH,
    show_progress: bool = True,
) -> None:
    """
    根据 DataFrame 拆分基因组序列.

    期望列: Chromosome, start, end, id
    """
    fasta = Fasta(fasta_path)
    iterator = bed_df.itertuples(index=False)
    total = len(bed_df)
    if show_progress:
        iterator = tqdm(iterator, total=total, desc="导出拆分后的序列")

    with out_fasta_path.open("w") as out:
        for idx, row in enumerate(iterator, start=1):
            seqid = str(row.Chromosome)
            start = int(row.start)
            end = int(row.end)
            name = str(row.id)

            if seqid not in fasta:
                raise ValueError(f"{seqid} not found in genome — skipped ({name})")

            seq = str(fasta[seqid][start:end])
            out.write(f">{name}\n")
            for chunk in textwrap.wrap(seq, width=line_length):
                out.write(chunk + "\n")

            logger.info("[{}/{}] Wrote {} ({} bp)", idx, total, out_fasta_path, len(seq))

    logger.success(
        "✅ Genome splitting completed! {} fragments written to {}",
        len(bed_df),
        out_fasta_path,
    )


def load_fai(genome_fai: Path) -> pd.DataFrame:
    """Load FASTA index information."""
    df = pd.read_table(
        genome_fai,
        header=None,
        names=["Chromosome", "chrom_size"],
        usecols=[0, 1],
        sep="\t",
    )
    if (df["chrom_size"] <= 0).any():
        raise ValueError("FAI 文件中的染色体大小必须为正数。")
    return df


def compute_split_windows(fai_df: pd.DataFrame, split_size: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Compute split start/end windows and separate chromosomes by size."""
    df = fai_df.copy()
    df["split_start"] = df["chrom_size"].apply(
        lambda size: max(0.3, 1 - split_size / size) * size
    )
    df["split_end"] = df["chrom_size"] - df["split_start"]
    need_split = df[df["chrom_size"] >= split_size].copy()
    no_split = df[df["chrom_size"] < split_size].copy()
    return need_split, no_split


def finalize_unsplit_chromosomes(df: pd.DataFrame) -> pd.DataFrame:
    """Prepare BED entries for chromosomes that do not require splitting."""
    if df.empty:
        return pd.DataFrame(columns=BED_COLUMNS)
    result = df.copy()
    result["start"] = 0
    result["end"] = result["chrom_size"]
    result["id"] = result["Chromosome"]
    return result[BED_COLUMNS]


def main(
    genome_fa: Path,
    genome_fai: Path,
    gff: Path,
    split_cat_bed: Path,
    min_gene_gap: int = 10_000,
    split_size: int = SPLIT_SIZE,
) -> None:
    fai_df = load_fai(genome_fai)
    need_to_split_df, do_not_need_to_split_df = compute_split_windows(
        fai_df, split_size
    )

    gene_df = parse_gff_with_pandas(gff)
    gap_df = calculate_gaps(gene_df)
    best_split_sites = select_split_candidates(need_to_split_df, gap_df, min_gene_gap)
    split_chr_bed = generate_split_chr_bed(best_split_sites)
    untouched_chr_bed = finalize_unsplit_chromosomes(do_not_need_to_split_df)
    split_chr_bed_all = pd.concat(
        [split_chr_bed, untouched_chr_bed], ignore_index=True
    )

    split_chr_bed_all.to_csv(
        split_cat_bed,
        sep="\t",
        index=False,
        header=False,
        columns=BED_COLUMNS,
    )

    out_gff = gff.with_suffix(".split.gff")
    logger.info("生成分割后的GFF文件: {}", out_gff)
    generate_split_gff(split_chr_bed, gff, out_gff)

    out_fasta = genome_fa.with_suffix(".split.fa")
    logger.info("生成分割后的基因组文件: {}", out_fasta)
    generate_split_genome(genome_fa, split_chr_bed_all, out_fasta)


if __name__ == "__main__":
    typer.run(main)
