from pathlib import Path

import delegator
import typer
from loguru import logger


def make_chain(
    ref: Path, query: Path, vcf: Path, outdir: Path, threads: int = 16
) -> None:
    """Make a chain file from a reference and query file."""
    ref_name = ref.stem
    query_name = query.stem
    outdir.mkdir(exist_ok=True)
    out_paf = outdir / f"{ref_name}_to_{query_name}.paf"
    minimap_cmd = f"minimap2 -t {threads} -cx asm5 --cs {query} {ref} > {out_paf}"
    logger.info(f"Running minimap2: {minimap_cmd}")
    delegator.run(minimap_cmd)
    chain_file = outdir / f"{ref_name}_to_{query_name}.chain"
    paf2chain_cmd = f"transanno minimap2chain {out_paf} --output {chain_file}"
    logger.info(f"Running paf2chain: {paf2chain_cmd}")
    delegator.run(paf2chain_cmd)
    lift_over_vcf = outdir / f"liftover.{ref_name}_to_{query_name}.vcf.gz"
    rejected_vcf = outdir / f"rejected.{ref_name}_to_{query_name}.vcf.gz"
    liftover_cmd = (
        f"transanno liftvcf --original-assembly {ref} "
        f"--new-assembly {query} "
        f"--chain {chain_file} --vcf {vcf} --output {lift_over_vcf} --fail {rejected_vcf} "
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
