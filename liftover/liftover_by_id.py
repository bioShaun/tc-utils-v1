from pathlib import Path

import delegator
import pandas as pd
import typer
from loguru import logger
from pyfaidx import Fasta
from tqdm import tqdm


def make_id_vcf(id_file: Path, ref_fa: Path, force: bool) -> Path:
    ref_fasta = Fasta(ref_fa)
    id_vcf_file = Path(f"{id_file}.vcf")
    if id_vcf_file.exists() and not force:
        logger.info(f"VCF file {id_vcf_file} already exists. Skipping creation.")
        return id_vcf_file
    with open(id_file, "r") as id_inf, open(id_vcf_file, "w") as vcf_inf:
        vcf_inf.write(f"##fileformat=VCFv4.2\n")
        vcf_inf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for line in tqdm(id_inf, desc="Processing IDs"):
            each_id = line.strip()
            chrom = "-".join(each_id.split("_")[:-1])
            pos = int(each_id.split("_")[-1])
            ref_seq = fetch_ref_nucleotide(ref_fasta, chrom, pos)
            vcf_inf.write(f"{chrom}\t{pos}\t{each_id}\t{ref_seq}\t.\t.\t.\t.\n")
    return id_vcf_file


def fetch_ref_nucleotide(ref_fa: Fasta, chrom: str, pos: int) -> str:
    return str(ref_fa[chrom][pos - 1 : pos].seq) if chrom in ref_fa else "N"


def make_chain(
    id_file: Path,
    chain: Path,
    ref_fa: Path,
    query_fa: Path,
    outdir: Path,
    force: bool = False,
) -> None:
    """Make a chain file from a reference and query file."""
    logger.info(f"Generating vcf from ID file: {id_file}")
    vcf: Path = make_id_vcf(id_file, ref_fa, force=force)
    lift_over_vcf: Path = outdir / f"liftover.{id_file.name}.vcf.gz"
    rejected_vcf: Path = outdir / f"rejected.{id_file.name}.vcf.gz"
    liftover_cmd: str = (
        f"transanno liftvcf --original-assembly {ref_fa} "
        f"--new-assembly {query_fa} "
        f"--chain {chain} --vcf {vcf} --output {lift_over_vcf} --fail {rejected_vcf} "
    )
    if force or not lift_over_vcf.is_file():
        logger.info(f"Running liftover command: {liftover_cmd}")
        delegator.run(liftover_cmd)
    else:
        logger.info(
            f"Lifted over VCF already exists: {lift_over_vcf}. Skipping liftover."
        )
    lift_bed = pd.read_table(
        lift_over_vcf,
        header=None,
        sep="\t",
        comment="#",
        usecols=[0, 1],
        names=["chrom", "pos"],
        compression="gzip",
    )
    lift_bed["start"] = lift_bed["pos"] - 1
    lift_bed["id"] = lift_bed["chrom"] + "_" + lift_bed["start"].astype(str)
    lift_bed.drop_duplicates(inplace=True)
    outdir.mkdir(exist_ok=True, parents=True)
    raw_bed = outdir / f"raw.{id_file.stem}.bed"
    probe_name = id_file.stem

    probe_id_file = outdir / f"{probe_name}.id"
    lift_bed.to_csv(probe_id_file, sep="\t", index=False, header=False, columns=["id"])
    lift_bed.to_csv(
        raw_bed, sep="\t", index=False, header=False, columns=["chrom", "start", "pos"]
    )
    ref_fa_idx = f"{ref_fa}.fai"
    probe_bed = outdir / f"{probe_name}.bed"
    sort_cmd = f"bedtools sort -i {raw_bed} -g {ref_fa_idx} > {probe_bed}"
    logger.info(f"Running bedtools sort: {sort_cmd}")
    delegator.run(sort_cmd)
    snp_calling_bed = outdir / f"{probe_name}.snpcalling.bed"
    span_bed_cmd = f"bedtools slop -i {probe_bed} -g {ref_fa_idx} -b 100 | bedtools merge -i - > {snp_calling_bed}"
    logger.info(f"Running bedtools span bed: {span_bed_cmd}")
    delegator.run(span_bed_cmd)
    # rm raw bed
    raw_bed.unlink()


if __name__ == "__main__":
    typer.run(make_chain)
