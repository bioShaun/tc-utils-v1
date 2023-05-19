from enum import Enum
import typer
import delegator
import pandas as pd
import numpy as np
from functools import partial, reduce
from operator import itemgetter
from loguru import logger

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

from pathlib import Path
from pandarallel import pandarallel

RAW_HEADER = ["chrom", "pos", "allele", "score"]

FLAT_HEADER = ["id", "chrom", "pos", "allele", "score"]

PROB_LENGTH = 120

MATCH_CUTOFFS = [0, 60, 100]


app = typer.Typer()


class InputType(str, Enum):
    TABLE = "table"
    VCF = "vcf"

    def __str__(self) -> str:
        return self.value


@app.command()
def flat(csv_input: Path) -> None:
    df = pd.read_csv(csv_input, header=None, names=RAW_HEADER)
    df.dropna(inplace=True)
    df["ref"] = df["allele"].map(lambda x: x.split(",")[0])
    df["alt"] = df["allele"].map(lambda x: x.split(",")[1:])
    df["alt_index"] = df["allele"].map(lambda x: list(range(1, len(x.split(",")))))
    df = df.explode(["alt", "alt_index"])
    df["id"] = df.apply(lambda x: f"{x['chrom']}_{x['pos']}_{x['alt_index']}", axis=1)

    df.drop("allele", axis=1, inplace=True)
    csv_out = csv_input.with_suffix(".flat.csv")
    df.to_csv(csv_out, index=False, columns=FLAT_HEADER)


def isIndel(allele: str) -> bool:
    if allele:
        seq_list = list(map(len, allele.split(",")))
        return max(seq_list) > 1
    return False


def get_seq(row, record, gc_bias):
    pos = row["pos"] - 1
    seq_list = []
    center_probe = record.seq[pos - PROB_LENGTH // 2 : pos + PROB_LENGTH // 2]
    center_probe_gc = gc_fraction(center_probe)
    probe_gc_bias = abs(center_probe_gc - 0.5)
    if center_probe.lower().count("n") > 10 or probe_gc_bias > gc_bias:
        center_candidates = []
        for i in range(10, PROB_LENGTH // 4, 10):
            tmp_probe1 = record.seq[
                pos - PROB_LENGTH // 2 - i : pos + PROB_LENGTH // 2 - i
            ]
            tmp_probe2 = record.seq[
                pos - PROB_LENGTH // 2 + i : pos + PROB_LENGTH // 2 + i
            ]
            tmp_probe1_n = tmp_probe1.lower().count("n")
            tmp_probe1_gc_bias = abs(gc_fraction(tmp_probe1) - 0.5)
            center_candidates.append([tmp_probe1, tmp_probe1_gc_bias, tmp_probe1_n])
            tmp_probe2_n = tmp_probe2.lower().count("n")
            tmp_probe2_gc_bias = abs(gc_fraction(tmp_probe2) - 0.5)
            center_candidates.append([tmp_probe2, tmp_probe2_gc_bias, tmp_probe2_n])
        center_candidates_sort = sorted(center_candidates, key=itemgetter(2, 1))
        seq_list.append(center_candidates_sort[0][0])
    else:
        seq_list.append(center_probe)

    if isIndel(row["allele"]):
        seq_list.append(record.seq[pos - PROB_LENGTH : pos])
        seq_list.append(record.seq[pos + 1 : pos + 1 + PROB_LENGTH])
    return seq_list


def get_seq_label(row):
    if isIndel(row["allele"]):
        return ["center", "left", "right"]
    return ["center"]


def get_ref(row, record):
    pos = row["pos"] - 1
    ref_len = len(row["ref"])
    seq = record.seq[pos : pos + ref_len]
    return str(seq)


def load_table_from_vcf(vcf: Path) -> pd.DataFrame:
    df = pd.read_csv(
        vcf,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 1, 3, 4],
        names=["chrom", "pos", "ref", "alt"],
    )
    df["allele"] = df.apply(lambda x: f"{x['ref']}>{x['alt']}", axis=1)
    df["score"] = 0
    return df


@app.command()
def addSeq(
    input_file: Path,
    ref: Path,
    input_type: InputType = InputType.TABLE,
    id_map: Path = typer.Option(None),
    gc_bias: float = 0.15,
    threads: int = 4,
) -> None:
    pandarallel.initialize(progress_bar=True, nb_workers=threads)
    if input_type == InputType.VCF:
        df = load_table_from_vcf(input_file)
    else:
        df = pd.read_csv(input_file, header=None, names=RAW_HEADER)
        df.allele.fillna("", inplace=True)
        dup_df = df[df.duplicated(subset=["chrom", "pos"])]
        if len(dup_df) > 0:
            logger.warning(f"Duplicate position found: {dup_df}")
            dup_df.to_csv(input_file.with_suffix(".dup.csv"), index=False)
            df.drop_duplicates(subset=["chrom", "pos"], inplace=True)
    id_map_dict = {}
    if id_map:
        id_map_dict = {line.split()[0]: line.split()[1] for line in id_map.open()}

    add_seq_df_list = []
    for record in SeqIO.parse(ref, "fasta"):
        chrom = record.id
        logger.info(f"add sequence for chromosome {chrom}")
        if id_map:
            if record.id not in id_map_dict:
                continue
            chrom = id_map_dict[record.id]
        chrom_df = df[df["chrom"] == chrom].copy()
        if chrom_df.empty:
            continue
        get_seq_by_chrom = partial(get_seq, record=record, gc_bias=gc_bias)
        chrom_df["sequence"] = list(chrom_df.parallel_apply(get_seq_by_chrom, axis=1))  # type: ignore
        chrom_df["sequence_type"] = list(chrom_df.parallel_apply(get_seq_label, axis=1))  # type: ignore
        add_seq_df_list.append(chrom_df)
    add_seq_df = pd.concat(add_seq_df_list)
    add_seq_df = add_seq_df.explode(["sequence", "sequence_type"])
    add_seq_df["GC"] = add_seq_df["sequence"].map(gc_fraction)
    add_seq_df["N_count"] = add_seq_df["sequence"].map(lambda x: x.lower().count("n"))
    add_seq_df["id"] = add_seq_df.parallel_apply(
        lambda x: f"{x['chrom']}_{x['pos']}_{x['sequence_type']}", axis=1
    )  # type: ignore
    csv_out = input_file.with_suffix(".seq.csv")
    add_seq_df.to_csv(
        csv_out,
        index=False,
        columns=FLAT_HEADER + ["sequence", "sequence_type", "GC", "N_count"],
        float_format="%.3f",
    )
    seq_file = input_file.with_suffix(".seq.fa")
    with seq_file.open("w") as f:
        for _, row in add_seq_df.iterrows():
            f.write(f">{row['id']}\n{row['sequence']}\n")


@app.command()
def add_seq_from_vcf(vcf: Path, ref: Path, id_map: Path = typer.Option(None)) -> None:
    df = pd.read_csv(
        vcf,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 1, 3, 4],
        names=["chrom", "pos", "ref", "alt"],
    )
    df["allele"] = df.apply(lambda x: f"{x['ref']}>{x['alt']}", axis=1)
    df["score"] = 0


@app.command()
def check_ref(csv_input: Path, ref: Path, id_map: Path = typer.Option(None)) -> None:
    df = pd.read_csv(csv_input)
    id_map_dict = {}
    if id_map:
        id_map_dict = {line.split()[0]: line.split()[1] for line in id_map.open()}

    notMatch_seq_df_list = []
    for record in SeqIO.parse(ref, "fasta"):
        chrom = record.id
        if id_map:
            if record.id not in id_map_dict:
                continue
            chrom = id_map_dict[record.id]
        chrom_df = df[df["chrom"] == chrom].copy()
        if chrom_df.empty:
            continue
        get_ref_by_chrom = partial(get_ref, record=record)
        chrom_df["check_ref"] = list(chrom_df.apply(get_ref_by_chrom, axis=1))
        notMatch_chrom_seq_df = chrom_df[chrom_df["ref"] != chrom_df["check_ref"]]
        notMatch_seq_df_list.append(notMatch_chrom_seq_df)
    notMatch_seq_df = pd.concat(notMatch_seq_df_list)
    csv_out = csv_input.with_suffix(".notMatch.csv")
    notMatch_seq_df.to_csv(csv_out, index=False, columns=FLAT_HEADER + ["check_ref"])


@app.command()
def addmaf(csv_input: Path, vcftools_frq: Path, threads: int = 4) -> None:
    pandarallel.initialize(progress_bar=True, nb_workers=threads)
    candidates_df = pd.read_csv(csv_input)
    frq_df = pd.read_csv(
        vcftools_frq,
        sep="\t",
        usecols=[0, 1, 4],
        skiprows=1,
        header=None,
        names=["chrom", "pos", "ref_freq"],
    )
    candidates_df = candidates_df.merge(frq_df)
    candidates_df["maf"] = candidates_df["ref_freq"].parallel_map(
        lambda x: x if x <= 0.5 else 1 - x
    )
    candidates_df.drop("ref_freq", inplace=True, axis=1)
    csv_out = csv_input.with_suffix(".maf.csv")
    candidates_df.to_csv(csv_out, index=False, float_format="%.3f")


@app.command()
def addMatch(
    flat_file: Path,
    seq_file: Path,
    ref_file: Path,
    minimap2: str = "minimap2",
    label: str = "candidates",
) -> None:
    paf_file = seq_file.with_suffix(".paf")
    if not paf_file.exists():
        logger.info(f"Run minimap2 to generate {paf_file}")
        delegator.run(f"{minimap2} -cx sr {ref_file} {seq_file} > {paf_file}")
    paf_df = pd.read_csv(paf_file, sep="\t", usecols=[0, 9], names=["id", "matches"])
    match_df = (
        paf_df.id.value_counts().reset_index().rename(columns={"count": "matches"})
    )
    flat_df = pd.read_csv(flat_file)
    merged_df_list = [flat_df]
    for cutoff in MATCH_CUTOFFS:
        cutoff_name = f"match_{cutoff}"
        filter_paf_df = paf_df[paf_df["matches"] >= cutoff]
        match_df = (
            filter_paf_df.id.value_counts()
            .reset_index()
            .rename(columns={"count": cutoff_name})
        )
        merged_df_list.append(match_df)
    merged_df = reduce(
        lambda left, right: pd.merge(left, right, on="id", how="outer"), merged_df_list
    )
    # merged_df = pd.merge(flat_df, match_df, on="id", how="outer")
    merged_df.fillna(0, inplace=True)
    merged_df["label"] = label

    csv_out = flat_file.with_suffix(".match.csv")
    merged_df.to_csv(csv_out, index=False)


@app.command()
def addLociMatch(
    flat_file: Path,
    seq_file: Path,
    ref_file: Path,
    minimap2: str = "minimap2",
    label: str = "candidates",
) -> None:
    paf_file = seq_file.with_suffix(".paf")
    if not paf_file.exists():
        logger.info(f"Run minimap2 to generate {paf_file}")
        delegator.run(f"{minimap2} -cx sr {ref_file} {seq_file} > {paf_file}")
    paf_df = pd.read_csv(
        paf_file,
        sep="\t",
        usecols=[0, 5, 7, 8, 9, 11],
        names=["id", "new_chrom", "new_start", "new_end", "matches", "mapq"],
    )
    new_loci_df = paf_df[
        (paf_df["mapq"] == 60) & (paf_df["matches"] >= PROB_LENGTH * 0.9)
    ].copy()
    new_loci_df["new_pos"] = new_loci_df["new_start"] + PROB_LENGTH // 2
    new_loci_df["new_id"] = new_loci_df.apply(
        lambda x: f"{x['new_chrom']}_{x['new_pos']}_center", axis=1
    )
    new_loci_df = new_loci_df[["id", "new_id", "new_chrom", "new_pos"]]
    id_map_df = new_loci_df[["id", "new_id"]]
    id_map_df.to_csv(flat_file.with_suffix(".id_map.csv"), index=False)

    match_df = (
        paf_df.id.value_counts().reset_index().rename(columns={"count": "matches"})
    )
    flat_df = pd.read_csv(flat_file)
    new_loci_df = new_loci_df.merge(flat_df)

    merged_df_list = [new_loci_df]
    for cutoff in MATCH_CUTOFFS:
        cutoff_name = f"match_{cutoff}"
        filter_paf_df = paf_df[paf_df["matches"] >= cutoff]
        match_df = (
            filter_paf_df.id.value_counts()
            .reset_index()
            .rename(columns={"count": cutoff_name})
        )
        merged_df_list.append(match_df)
    merged_df = reduce(
        lambda left, right: pd.merge(left, right, on="id", how="left"), merged_df_list
    )
    # merged_df = pd.merge(flat_df, match_df, on="id", how="outer")
    merged_df.drop(["id", "chrom", "pos"], axis=1, inplace=True)
    merged_df.rename(
        columns={"new_chrom": "chrom", "new_pos": "pos", "new_id": "id"}, inplace=True
    )
    merged_df.fillna(0, inplace=True)
    merged_df["label"] = label

    csv_out = flat_file.with_suffix(".match.csv")
    merged_df.to_csv(csv_out, index=False)


if __name__ == "__main__":
    app()
