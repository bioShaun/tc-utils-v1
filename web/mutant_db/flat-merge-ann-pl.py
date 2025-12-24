from pathlib import Path
from typing import Tuple

import pandas as pd
import typer
import polars as pl


def flat_gt_table(gt_table: Path) -> Path:
    flat_gt_file = gt_table.with_suffix(".flat.tsv")
    with (
        open(flat_gt_file, "w") as flat_inf,
        open(gt_table) as gt_inf
        ):
        for line in gt_inf:
            chrom, pos, refer, alt, *score_and_sample = line.strip().split("\t")
            try:
                for sample_i in score_and_sample[1:]:
                    flat_inf.write(f"{chrom}\t{pos}\t{refer}\t{alt}\t{sample_i}\n")
            except ValueError:
                print(line)
                raise ValueError("sample list wrong type")
    return flat_gt_file


def split_accession(accession: str) -> Tuple[str, str, int, int]:
    sample_id, allele_info = accession.split("=")
    genotype, allele_depth = allele_info.split(";")
    if "." in genotype:
        genotype = "./."
        reference_depth, alternate_depth = 0, 0
    else:
        allele_depth = allele_depth.replace(".", "0")
        reference_depth, alternate_depth = allele_depth.split(",")
    return sample_id, genotype, int(reference_depth), int(alternate_depth)


def main(gt_table: Path, ann_table: Path, out_table: Path):
    flat_gt_table_file = flat_gt_table(gt_table)
    flat_gt_df = pl.scan_csv(
        flat_gt_table_file, 
        separator="\t", 
        has_header=False,
        new_columns=["chrom", "pos", "refer", "alt", "accession"]
    ).with_columns(
        pl.col("chrom").cast(pl.Utf8),
        pl.col("pos").cast(pl.Int64),
        pl.col("refer").cast(pl.Utf8),
        pl.col("alt").cast(pl.Utf8),
        pl.col("accession").cast(pl.Utf8),
    ).filter(~pl.col("alt").str.contains(","))

    flat_gt_df = flat_gt_df.with_columns(
        variant=pl.col("chrom").str.cat(pl.col("pos").cast(pl.Utf8), separator="_"),
        sample_id=pl.col("accession").str.split("=").list.get(0),
        genotype=pl.col("accession").str.split("=").list.get(1).str.split(";").list.get(0),
        ref_depth=pl.col("accession").str.split("=").list.get(1).str.split(";").list.get(1).str.split(",").list.get(0).cast(pl.Int64),
        alt_depth=pl.col("accession").str.split("=").list.get(1).str.split(";").list.get(1).str.split(",").list.get(1).cast(pl.Int64),
    ).filter(pl.col("genotype") != "./.")
    print(flat_gt_df)
    ann_df = pl.scan_csv(
        ann_table,
        separator="\t",
        has_header=False,
        new_columns=[
            "chrom",
            "pos",
            "refer",
            "alt",
            "type",
            "impact",
            "gene",
            "exon_rank",
            "cds_pos",
            "protein_pos",
        ],
    ).select(pl.all().head(10))
    print(ann_df)
    ann_df = ann_df.with_columns(
        chrom=pl.col("chrom").cast(pl.Utf8),
        pos=pl.col("pos").cast(pl.Int64),
        refer=pl.col("refer").cast(pl.Utf8),
        alt=pl.col("alt").cast(pl.Utf8),
    ).filter(~pl.col("alt").str.contains(","))
    ann_df = ann_df.unique(subset=["chrom", "pos", "refer", "alt"])
    add_ann_df = ann_df.join(flat_gt_df, on=["chrom", "pos", "refer", "alt"]).drop("accession")
    add_ann_df.collect().write_csv(out_table, separator="\t")


if __name__ == "__main__":
    typer.run(main)
