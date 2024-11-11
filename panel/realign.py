import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import List, Optional, Tuple

import delegator
import numpy as np
import pandas as pd
import typer
from loguru import logger
from tqdm import tqdm

FLANK_SIZE = 60


app = typer.Typer()


class TargetType(str, Enum):
    bed = "bed"
    vcf = "vcf"


@dataclass
class InsDelCount:
    ins_count: int
    del_count: int


def get_cigar_list(cigar: str) -> List[Tuple[int, str]]:
    """
    Parse a CIGAR string into a list of tuples.

    Args:
    cigar (str): A CIGAR string (e.g., "3M1I4M2D5M")

    Returns:
    List[Tuple[int, str]]: A list of tuples where each tuple contains
                           (count, operation)

    Example:
    >>> get_cigar_list("3M1I4M2D5M")
    [(3, 'M'), (1, 'I'), (4, 'M'), (2, 'D'), (5, 'M')]
    """
    # Use a single regex to split the CIGAR string into pairs of (count, operation)
    cigar_pairs = re.findall(r"(\d+)([MIDNSHPX=])", cigar)

    # Convert the count to integer and return the list of tuples
    return [(int(count), operation) for count, operation in cigar_pairs]


def get_del_ins(row: pd.Series) -> Optional[InsDelCount]:
    """
    Analyze CIGAR string to count insertions and deletions up to a certain offset.

    Args:
    row (pd.Series): A pandas Series containing 'cigar', 'offset_start', 'offset_end', and 'strand' information.

    Returns:
    Optional[InsDelCount]: A named tuple with insertion and deletion counts, or None if the offset is not reached.

    Raises:
    ValueError: If an unexpected CIGAR operation is encountered.
    """
    cigar_info = get_cigar_list(row["cigar"])
    offset_from_target = (
        row["offset_start"] if row["strand"] == "+" else row["offset_end"]
    )

    match_count = ins_count = del_count = 0

    for count, operation in cigar_info:
        if operation == "M":
            match_count += count
        elif operation == "I":
            ins_count += count
        elif operation == "D":
            if match_count == offset_from_target:
                return None
            del_count += count
        else:
            raise ValueError(f"Unexpected CIGAR operation: {operation}")

        if match_count > offset_from_target:
            break

    if match_count < offset_from_target:
        return None  # Offset not reached

    return InsDelCount(ins_count=ins_count, del_count=del_count)


def get_pos(row: pd.Series) -> Optional[int]:
    indel_ins_info = get_del_ins(row)
    if indel_ins_info is None:
        return None
    indel_bias = indel_ins_info.del_count - indel_ins_info.ins_count
    if row["strand"] == "+":
        return (
            row["match_start"]
            + row["offset_start"]
            + 1
            - row["probe_start"]
            + indel_bias
        )
    return row["match_start"] + row["offset_end"] + 1 - row["probe_start"] + indel_bias


def paf2idmap(
    paf_df: pd.DataFrame,
    match_cutoff: float,
    offset_df: pd.DataFrame,
    out_prefix: Path,
    save_file: bool = True,
) -> Optional[pd.DataFrame]:

    paf_df.sort_values(
        [
            "match_length",
            "mapq",
        ],
        ascending=False,
        inplace=True,
    )
    paf_df.drop_duplicates(subset=["id"], inplace=True)
    paf_df["match_ratio"] = paf_df["match_length"] / paf_df["probe_length"]
    filter_df = paf_df[paf_df["match_ratio"] > match_cutoff].copy()
    filter_df = filter_df.merge(offset_df, on="id")
    filter_df["pos"] = filter_df.apply(get_pos, axis=1)
    filter_df.dropna(subset=["pos"], inplace=True)
    filter_df["pos"] = filter_df["pos"].astype("int")
    filter_df["new_id"] = filter_df.apply(lambda x: f'{x["chrom"]}_{x["pos"]}', axis=1)
    filter_df["pos_0"] = filter_df["pos"] - 1
    id_map = out_prefix.with_suffix(".idmap.tsv")
    id_map_df = filter_df[["id", "new_id"]].copy()
    if save_file:
        filter_df.to_csv(
            id_map, header=False, index=False, columns=["id", "new_id"], sep="\t"
        )
        id_map_target_bed = id_map.with_suffix(".target.bed")
        filter_df.sort_values(["chrom", "pos"], inplace=True)
        filter_df.to_csv(
            id_map_target_bed,
            sep="\t",
            index=False,
            header=False,
            columns=["chrom", "pos_0", "pos", "id"],
        )

    return id_map_df


def generate_bed_from_vcf(vcf: Path) -> Path:
    vcf_df = pd.read_table(
        vcf,
        header=None,
        usecols=[
            0,
            1,
        ],
        comment="#",
        names=["chrom", "end"],
    )
    vcf_df["start"] = vcf_df["end"] - 1
    vcf_df["id"] = vcf_df.apply(lambda x: f'{x["chrom"]}_{x["end"]}', axis=1)
    bed_file = vcf.with_suffix(".bed")
    vcf_df.to_csv(
        bed_file,
        sep="\t",
        index=False,
        header=False,
        columns=["chrom", "start", "end", "id"],
    )
    return bed_file


def generate_flank_bed(target_bed: Path, flank_size: int, genome_fai: Path) -> Path:
    flank_bed = target_bed.with_suffix(f".flank{flank_size}.bed")
    cmd_line = (
        f"bedtools slop -b {flank_size} -i {target_bed} -g {genome_fai} > {flank_bed}"
    )
    delegator.run(cmd_line)
    return flank_bed


def generate_flank_fa(flank_bed: Path, genome_fa: Path, force: bool) -> Path:
    flank_fa = flank_bed.with_suffix(".fa")
    if force or not flank_fa.is_file():
        cmd_line = f"bedtools getfasta -fi {genome_fa} -fo {flank_fa} -bed {flank_bed} -nameOnly"
        delegator.run(cmd_line)
    return flank_fa


def generate_flank_paf(
    flank_fa: Path, genome_sr_idx: Path, threads: int, force: bool
) -> Path:
    flank_paf = flank_fa.with_suffix(".paf")
    if force or not flank_paf.is_file():
        cmd_line = (
            f"minimap2 -t {threads} -cx sr {genome_sr_idx} {flank_fa} > {flank_paf}"
        )
        delegator.run(cmd_line)
    return flank_paf


def generate_offset_df(target_bed: Path, flank_bed: Path) -> pd.DataFrame:
    target_df = pd.read_table(
        target_bed, header=None, names=["target_start", "id"], usecols=[1, 3]
    )
    flank_df = pd.read_table(
        flank_bed,
        header=None,
        names=["flank_start", "flank_end", "id"],
        usecols=[1, 2, 3],
    )
    offset_df = target_df.merge(flank_df)
    offset_df["offset_start"] = offset_df["target_start"] - offset_df["flank_start"]
    offset_df["offset_end"] = offset_df["flank_end"] - offset_df["target_start"] - 1
    offset_df.drop(columns=["target_start", "flank_start", "flank_end"], inplace=True)
    return offset_df


def probe_sequence_from_flank(flank: str) -> str:
    left = flank.split("[")[0]
    center = flank.split("[")[1][0]
    right = flank.split("]")[1]
    return f"{left}{center}{right}"


def fasta_from_probe_table(probe_table: Path) -> Tuple[pd.DataFrame, Path]:
    probe_df = pd.read_table(probe_table)
    probe_fasta = probe_table.with_suffix(".fa")
    probe_df["sequence"] = probe_df["Flank"].map(probe_sequence_from_flank)
    with open(probe_fasta, "w", encoding="utf-8") as f:
        for row in probe_df.itertuples():
            f.write(f">{row.id}\n{row.sequence}\n")
    probe_df["seq_length"] = probe_df["sequence"].map(len)
    probe_df["offset_start"] = probe_df["Flank"].map(lambda x: x.index("["))
    probe_df["offset_end"] = probe_df["seq_length"] - probe_df["offset_start"] - 1
    return probe_df[["id", "offset_start", "offset_end"]].copy(), probe_fasta


def extract_cigar_from_line(line: str) -> str:
    """Extract the cigar string from a Minimap2 alignment line."""
    match = re.search(r"cg:Z:(\w+)", line)
    if match:
        return match.group(1)
    raise ValueError(f"cigar not found in line: {line!r}")


def cigar_list_from_paf(paf: Path) -> pd.DataFrame:
    cigar_list = [extract_cigar_from_line(each) for each in paf.open()]
    return pd.DataFrame(cigar_list, columns=["cigar"])


@app.command()
def realign(
    target_file: Path,
    genome: Path,
    genome_sr_idx: Path,
    threads: int = 16,
    cut_off: float = 0.5,
    force: bool = False,
    target_type: TargetType = TargetType.bed,
) -> None:
    """
    Main function to perform a series of operations based on the input parameters.

    Parameters:
        target_bed (Path): Path to the target bed file.
        genome (Path): Path to the genome file.
        genome_sr_idx (Path): Path to the genome index file.
        id_map (Path): Path to the ID mapping file.
        threads (int, optional): Number of threads to use. Defaults to 16.

    Returns:
        None
    """
    genome_fai = genome.parent / f"{genome.name}.fai"
    if target_type == TargetType.vcf:
        logger.info(f"Generating bed file from vcf...")
        target_bed = generate_bed_from_vcf(target_file)
    else:
        target_bed = target_file
    logger.info(f"Generating {FLANK_SIZE} bp flanks bed...")
    flank_bed = generate_flank_bed(
        target_bed=target_bed, flank_size=FLANK_SIZE, genome_fai=genome_fai
    )
    logger.info(f"Generating {FLANK_SIZE} bp flanks fasta...")
    flank_fa = generate_flank_fa(
        flank_bed=flank_bed,
        genome_fa=genome,
        force=force,
    )
    logger.info(f"Generating {FLANK_SIZE} bp flanks paf...")
    flank_paf = generate_flank_paf(
        flank_fa=flank_fa, genome_sr_idx=genome_sr_idx, threads=threads, force=force
    )
    offset_df = generate_offset_df(target_bed=target_bed, flank_bed=flank_bed)
    logger.info(f"Generating {FLANK_SIZE} bp flanks id map...")
    # match_cutoff = np.ceil(cut_off * 2 * FLANK_SIZE)
    paf_df = pd.read_table(
        flank_paf,
        header=None,
        usecols=[0, 1, 2, 4, 5, 7, 9, 11],
        names=[
            "id",
            "probe_length",
            "probe_start",
            "strand",
            "chrom",
            "match_start",
            "match_length",
            "mapq",
        ],
    )
    cigar_df = cigar_list_from_paf(flank_paf)
    add_cigar_paf_df = pd.concat([paf_df, cigar_df], axis=1)
    paf2idmap(
        paf_df=add_cigar_paf_df,
        offset_df=offset_df,
        match_cutoff=cut_off,
        out_prefix=flank_paf,
    )


@app.command()
def realign2(
    probe_table: Path,
    genome_sr_idx: Path,
    threads: int = 16,
    cut_off: float = 0.5,
    force: bool = False,
    allow_ins_del: bool = False,
) -> None:
    """
    Main function to perform a series of operations based on the input parameters.

    Parameters:
        fa (Path): Path to the fasta file.
        offset_bed (Path): Path to the offset table file.
        genome (Path): Path to the genome file.
        genome_sr_idx (Path): Path to the genome index file.
        id_map (Path): Path to the ID mapping file.
        threads (int, optional): Number of threads to use. Defaults to 16.

    Returns:
        None
    """
    logger.info(f"generate fa and offset ...")
    offset_df, flank_fa = fasta_from_probe_table(probe_table=probe_table)
    logger.info(f"map to genome ...")
    flank_paf = generate_flank_paf(
        flank_fa=flank_fa, genome_sr_idx=genome_sr_idx, threads=threads, force=force
    )
    logger.info(f"Generating id map...")
    paf_df = pd.read_table(
        flank_paf,
        header=None,
        usecols=[0, 1, 2, 4, 5, 7, 9, 11],
        names=[
            "id",
            "probe_length",
            "probe_start",
            "strand",
            "chrom",
            "match_start",
            "match_length",
            "mapq",
        ],
    )
    cigar_df = cigar_list_from_paf(flank_paf)
    add_cigar_paf_df = pd.concat([paf_df, cigar_df], axis=1)
    paf2idmap(
        paf_df=add_cigar_paf_df,
        offset_df=offset_df,
        match_cutoff=cut_off,
        out_prefix=flank_paf,
    )


@app.command()
def adjust_annotation_by_realign(annotation: Path, id_map: Path) -> None:
    anno_df = pd.read_table(annotation)
    id_map_df = pd.read_table(id_map, header=None, names=["id", "new_id"])
    re_id_df = id_map_df.merge(anno_df)
    re_id_df["chrom"] = re_id_df["new_id"].map(lambda x: "_".join(x.split("_")[:-1]))
    re_id_df["pos"] = re_id_df["new_id"].map(lambda x: int(x.split("_")[-1]))
    re_id_df.drop(columns=["new_id"], inplace=True)
    if "probe_start" in re_id_df.columns:
        re_id_df["probe_length"] = re_id_df["probe_end"] - re_id_df["probe_start"]
        re_id_df["probe_start"] = (
            re_id_df["pos"] - 1 - (re_id_df["probe_length"] // 2) + re_id_df["offset"]
        )
        re_id_df["probe_end"] = re_id_df["probe_start"] + re_id_df["probe_length"]
        re_id_df.drop(columns=["probe_length"], inplace=True)
    realign_annotation = annotation.with_suffix(".realign.tsv")
    re_id_df.to_csv(realign_annotation, sep="\t", index=False)


def blast2paf(blast_df: pd.DataFrame, probe_length: int) -> pd.DataFrame:
    paf_df = blast_df.copy()
    paf_df["probe_length"] = probe_length
    pos_strand_df = paf_df[paf_df["sstart"] < paf_df["send"]].copy()
    neg_strand_df = paf_df[paf_df["sstart"] > paf_df["send"]].copy()
    neg_strand_df["sstart"], neg_strand_df["send"] = (
        neg_strand_df["send"],
        neg_strand_df["sstart"],
    )
    pos_strand_df["qstart"] = pos_strand_df["qstart"] - 1
    neg_strand_df["qstart"] = probe_length - neg_strand_df["qend"]
    pos_strand_df["strand"] = "+"
    neg_strand_df["strand"] = "-"
    paf_df = pd.concat([pos_strand_df, neg_strand_df])
    paf_df = paf_df[
        [
            "qseqid",
            "probe_length",
            "qstart",
            "strand",
            "sseqid",
            "sstart",
            "qlength",
            "bitscore",
        ]
    ].copy()
    paf_df["sstart"] = paf_df["sstart"] - 1
    paf_df.rename(
        columns={
            "qseqid": "id",
            "qstart": "probe_start",
            "sseqid": "chrom",
            "sstart": "match_start",
            "qlength": "match_length",
            "bitscore": "mapq",
        },
        inplace=True,
    )
    return paf_df


@app.command()
def realign3(
    annotation: Path,
    blast_dir: Path,
    probe_length: int,
    miss_cut_off_bp: int,
) -> None:
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    # blast_dfs = [
    #     pd.read_table(
    #         each,
    #         header=None,
    #         names=[
    #             "qseqid",
    #             "sseqid",
    #             "pident",
    #             "qlength",
    #             "mismatch",
    #             "gapopen",
    #             "qstart",
    #             "qend",
    #             "sstart",
    #             "send",
    #             "evalue",
    #             "bitscore",
    #         ],
    #     )
    #     for each in blast_dir.glob("*")
    # ]

    blast_files = list(blast_dir.glob("*"))
    for n, blast_file_i in enumerate(tqdm(blast_files)):
        blast_df = pd.read_table(
            blast_file_i,
            header=None,
            names=[
                "qseqid",
                "sseqid",
                "pident",
                "qlength",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore",
            ],
        )
        match_cut_off_bp = probe_length - miss_cut_off_bp
        blast_df = blast_df[blast_df["qlength"] >= match_cut_off_bp].copy()
        blast_df = blast_df[blast_df["mismatch"] <= miss_cut_off_bp].copy()
        blast_df = blast_df[blast_df["gapopen"] == 0].copy()
        blast_df.drop_duplicates(subset=["qseqid"], inplace=True)
        paf_df = blast2paf(blast_df, probe_length)
        anno_df = pd.read_table(annotation)
        anno_df["offset_start"] = anno_df["pos"] - 1 - anno_df["probe_start"]
        anno_df["offset_end"] = probe_length - anno_df["offset_start"] - 1
        offset_table = anno_df[["id", "offset_start", "offset_end"]].copy()
        id_map_df = paf2idmap(
            paf_df=paf_df,
            match_cutoff=0,
            offset_df=offset_table,
            out_prefix=annotation,
            save_file=False,
        )
        mode = "w" if n == 0 else "a"
        if id_map_df is not None:
            id_map_df.to_csv(
                annotation.with_suffix(".idmap.tsv"),
                sep="\t",
                header=False,
                index=False,
                mode=mode,
            )


if __name__ == "__main__":
    app()
