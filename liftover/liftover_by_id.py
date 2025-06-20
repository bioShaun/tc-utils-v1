from pathlib import Path

import delegator
import pandas as pd
import typer
from loguru import logger
from pyfaidx import Fasta


def make_id_vcf(id_file: Path, ref_fa: Path) -> Path:
    ref_fasta = Fasta(ref_fa)
    id_vcf_file = Path(f"{id_file}.vcf")
    with open(id_file, "r") as id_inf, open(id_vcf_file, "w") as vcf_inf:
        vcf_inf.write(f"##fileformat=VCFv4.2\n")
        vcf_inf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for line in id_inf:
            each_id = line.strip()
            chrom = "-".join(each_id.split("_")[:-1])
            pos = each_id.split("_")[-1]
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
    out_dir: Path,
) -> None:
    """Make a chain file from a reference and query file."""
    vcf: Path = make_id_vcf(id_file)
    outdir: Path = id_file.parent
    lift_over_vcf: Path = outdir / f"liftover.{id_file.name}.vcf.gz"
    rejected_vcf: Path = outdir / f"rejected.{id_file.name}.vcf.gz"
    liftover_cmd: str = (
        f"transanno liftvcf --original-assembly {ref_fa} "
        f"--new-assembly {query_fa} "
        f"--chain {chain} --vcf {vcf} --output {lift_over_vcf} --fail {rejected_vcf} "
    )
    logger.info(f"Running liftover: {liftover_cmd}")
    delegator.run(liftover_cmd)
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
    out_dir.mkdir(exist_ok=True, parents=True)
    raw_bed = out_dir / f"raw.{id_file.stem}.bed"
    probe_name = id_file.stem

    probe_id_file = out_dir / f"{probe_name}.id"
    lift_bed.to_csv(probe_id_file, sep="\t", index=False, header=False, columns=["id"])
    lift_bed.to_csv(
        raw_bed, sep="\t", index=False, header=False, names=["chrom", "start", "pos"]
    )
    ref_fa_idx = f"{ref_fa}.fai"
    probe_bed = out_dir / f"{probe_name}.bed"
    sort_cmd = f"bedtools sort -i {raw_bed} -g {ref_fa_idx} > {probe_bed}"
    logger.info(f"Running bedtools sort: {sort_cmd}")
    delegator.run(sort_cmd)
    snp_calling_bed = out_dir / f"{probe_name}.snpcalling.bed"
    spa missed_sample_set = set() missed_sample_set = set()n_bed_cmd = f"bedtools slop -i {probe_bed} -g {ref_fa_idx} -b 200 | bedtools merge -i - > {snp_calling_bed}"
    logger.info(f"Running bedtools span bed: {span_bed_cmd}")
    delegator.run(span_bed_cmd)
    # rm raw bed
    raw_bed.unlink()


if __name__ == "__main__":
    typer.run(make_chain)


def assign_groups(node, threshold, group_dict=None):
    if group_dict is None:
        group_dict = {}

    if node.is_leaf():
        return group_dict

    dist = node.get_distance(node.get_tree_root())
    if dist > threshold:
        group_num = int(dist // threshold) + 1
        for leaf in node.get_leaves():
            # 只有当样品还没有被分组时才分配新组
            if leaf.name not in group_dict:
                group_dict[leaf.name] = f"Group{group_num}"

    for child in node.children:
        assign_groups(child, threshold, group_dict)

    return group_dict
