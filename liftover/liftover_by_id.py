from pathlib import Path

import delegator
import typer
from loguru import logger


def make_id_vcf(id_file: Path) -> Path:
    id_vcf_file = Path(f"{id_file}.vcf")
    with open(id_file, "r") as id_inf, open(id_vcf_file, "w") as vcf_inf:
        vcf_inf.write(f"##fileformat=VCFv4.2\n")
        for line in id_inf:
            each_id = line.strip()
            chrom = "-".join(each_id.split("_")[:-1])
            pos = each_id.split("_")[-1]
            vcf_inf.write(f"{chrom}\t{pos}\t{each_id}\t.\t.\t.\t.\t.\n")
    return id_vcf_file


def make_chain(
    id_file: Path,
    chain: Path,
    ref_fa: Path,
    query_fa: Path,
) -> None:
    """Make a chain file from a reference and query file."""
    vcf = make_id_vcf(id_file)
    outdir = id_file.parent
    lift_over_vcf = outdir / f"liftover.{id_file.name}.vcf.gz"
    rejected_vcf = outdir / f"rejected.{id_file.name}.vcf.gz"
    liftover_cmd = (
        f"transanno liftvcf --original-assembly {ref_fa} "
        f"--new-assembly {query_fa} "
        f"--chain {chain} --vcf {vcf} --output {lift_over_vcf} --fail {rejected_vcf} "
    )
    logger.info(f"Running liftover: {liftover_cmd}")
    delegator.run(liftover_cmd)


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
