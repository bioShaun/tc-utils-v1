from pathlib import Path
from typing import Tuple

import pandas as pd
import typer


def flat_gt_table(gt_table: Path) -> Path:
    flat_gt_file = gt_table.with_suffix(".flat.tsv")
    with open(flat_gt_file, "w") as flat_inf:
        with open(gt_table) as gt_inf:
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
        reference_depth, alternate_depth = allele_depth.split(",")
    return sample_id, genotype, int(reference_depth), int(alternate_depth)


def main(gt_table: Path, ann_table: Path, out_table: Path):
    flat_gt_table_file = flat_gt_table(gt_table)
    flat_gt_df = pd.read_table(
        flat_gt_table_file,
        header=None,
        names=["chrom", "pos", "refer", "alt", "accession"],
    )
    flat_gt_df = flat_gt_df[~flat_gt_df["alt"].str.contains(",")]
    flat_gt_df["variant"] = flat_gt_df["chrom"].str.cat(
        flat_gt_df["pos"].astype("str"), sep="_"
    )
    flat_gt_df[["sample_id", "genotype", "ref_depth", "alt_depth"]] = (
        flat_gt_df["accession"].map(split_accession).to_list()  # type: ignore
    )
    flat_gt_df = flat_gt_df[flat_gt_df["genotype"] != "./."]
    ann_df = pd.read_table(
        ann_table,
        header=None,
        names=[
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
        usecols=list(range(10)),
    )
    ann_df.drop_duplicates(subset=["chrom", "pos", "refer", "alt"], inplace=True)
    add_ann_df = ann_df.merge(flat_gt_df).drop("accession", axis=1)
    add_ann_df.to_csv(out_table, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)
