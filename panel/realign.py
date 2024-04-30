from pathlib import Path

import delegator
import numpy as np
import pandas as pd
import typer
from loguru import logger

FLANK_SIZE = 60


app = typer.Typer()


def get_pos(row: pd.Series) -> int:
    if row["strand"] == "+":
        return row["match_start"] + row["offset_start"] + 1 - row["probe_start"]
    return row["match_start"] + row["offset_end"] + 1 - row["probe_start"]


def paf2idmap(paf: Path, match_cutoff: int, offset_df: pd.DataFrame) -> None:
    paf_df = pd.read_table(
        paf,
        header=None,
        usecols=[0, 2, 4, 5, 7, 9, 11],
        names=[
            "id",
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


def generate_flank_bed(target_bed: Path, flank_size: int, genome_fai: Path) -> Path:
    flank_bed = target_bed.with_suffix(f".flank{flank_size}.bed")
    cmd_line = (
        f"bedtools slop -b {flank_size} -i {target_bed} -g {genome_fai} > {flank_bed}"
    )
    delegator.run(cmd_line)
    return flank_bed


def generate_flank_fa(flank_bed: Path, genome_fa: Path) -> Path:
    flank_fa = flank_bed.with_suffix(".fa")
    cmd_line = (
        f"bedtools getfasta -fi {genome_fa} -fo {flank_fa} -bed {flank_bed} -nameOnly"
    )
    delegator.run(cmd_line)
    return flank_fa


def generate_flank_paf(flank_fa: Path, genome_sr_idx: Path, threads: int) -> Path:
    flank_paf = flank_fa.with_suffix(".paf")
    cmd_line = f"minimap2 -t {threads} -cx sr {genome_sr_idx} {flank_fa} > {flank_paf}"
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


@app.command()
def realign(
    target_bed: Path,
    genome: Path,
    genome_sr_idx: Path,
    threads: int = 16,
    cut_off: float = 0.5,
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
    logger.info(f"Generating {FLANK_SIZE} bp flanks bed...")
    flank_bed = generate_flank_bed(
        target_bed=target_bed, flank_size=FLANK_SIZE, genome_fai=genome_fai
    )
    logger.info(f"Generating {FLANK_SIZE} bp flanks fasta...")
    flank_fa = generate_flank_fa(flank_bed=flank_bed, genome_fa=genome)
    logger.info(f"Generating {FLANK_SIZE} bp flanks paf...")
    flank_paf = generate_flank_paf(
        flank_fa=flank_fa, genome_sr_idx=genome_sr_idx, threads=threads
    )
    offset_df = generate_offset_df(target_bed=target_bed, flank_bed=flank_bed)
    logger.info(f"Generating {FLANK_SIZE} bp flanks id map...")
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
