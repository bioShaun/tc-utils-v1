from enum import Enum
from pathlib import Path
from typing import Tuple

import delegator
import numpy as np
import pandas as pd
import typer
from loguru import logger

FLANK_SIZE = 60


app = typer.Typer()


class TargetType(str, Enum):
    bed = "bed"
    vcf = "vcf"


def get_pos(row: pd.Series) -> int:
    if row["strand"] == "+":
        return row["match_start"] + row["offset_start"] + 1 - row["probe_start"]
    return row["match_start"] + row["offset_end"] + 1 - row["probe_start"]


def paf2idmap(paf: Path, match_cutoff: float, offset_df: pd.DataFrame) -> None:
    paf_df = pd.read_table(
        paf,
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
    filter_df = paf_df[paf_df["match_length"] > match_cutoff].copy()
    filter_df = filter_df.merge(offset_df, on="id")
    filter_df["pos"] = filter_df.apply(get_pos, axis=1)
    filter_df["new_id"] = filter_df.apply(lambda x: f'{x["chrom"]}_{x["pos"]}', axis=1)
    filter_df["pos_0"] = filter_df["pos"] - 1
    id_map = paf.with_suffix(".idmap.tsv")
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
        columns=["chrom", "pos_0", "pos"],
    )


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
    paf2idmap(paf=flank_paf, offset_df=offset_df, match_cutoff=cut_off)


@app.command()
def realign2(
    probe_table: Path,
    genome_sr_idx: Path,
    threads: int = 16,
    cut_off: float = 0.5,
    force: bool = False,
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
    match_cutoff = np.ceil(cut_off * 2 * FLANK_SIZE)
    paf2idmap(paf=flank_paf, offset_df=offset_df, match_cutoff=match_cutoff)


@app.command()
def adjust_annotation_by_realign(annotation: Path, id_map: Path) -> None:
    anno_df = pd.read_table(annotation)
    id_map_df = pd.read_table(id_map, header=None, names=["id", "new_id"])
    re_id_df = (
        id_map_df.merge(anno_df).drop("id", axis=1).rename(columns={"new_id": "id"})
    )
    re_id_df["chrom"] = re_id_df["id"].map(lambda x: "_".join(x.split("_")[:-1]))
    re_id_df["pos"] = re_id_df["id"].map(lambda x: int(x.split("_")[-1]))
    re_id_df["probe_length"] = re_id_df["probe_end"] - re_id_df["probe_start"]
    re_id_df["probe_start"] = (
        re_id_df["pos"] - 1 - (re_id_df["probe_length"] // 2) + re_id_df["offset"]
    )
    re_id_df["probe_end"] = re_id_df["probe_start"] + re_id_df["probe_length"]
    re_id_df.drop(columns=["probe_length"], inplace=True)
    realign_annotation = annotation.with_suffix(".realign.tsv")
    re_id_df.to_csv(realign_annotation, sep="\t", index=False)


if __name__ == "__main__":
    app()
