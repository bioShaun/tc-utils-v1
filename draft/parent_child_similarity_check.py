#!/usr/bin/env python3
import argparse
import pandas as pd
from cyvcf2 import VCF
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path


def parse_list(file_path):
    """读取样品列表"""
    if not file_path:
        return []
    return [x.strip() for x in open(file_path) if x.strip()]


def calc_similarity(vcf_path, child, parent):
    """计算单个子代-亲本配对的纯合位点相似度"""
    vcf = VCF(vcf_path, gts012=True)
    samples = vcf.samples

    try:
        idx_child = samples.index(child)
        idx_parent = samples.index(parent)
    except ValueError:
        return (child, parent, 0, 0, 0.0)

    n_sites = 0
    n_shared_hom = 0

    for variant in vcf:
        gt_child = variant.genotypes[idx_child][:2]
        gt_parent = variant.genotypes[idx_parent][:2]
        if -1 in gt_child or -1 in gt_parent:
            continue

        # 仅在两者均纯合时计算
        if gt_child[0] == gt_child[1] and gt_parent[0] == gt_parent[1]:
            n_sites += 1
            if gt_child[0] == gt_parent[0]:
                n_shared_hom += 1

    vcf.close()
    sim = n_shared_hom / n_sites if n_sites > 0 else 0
    return (child, parent, n_sites, n_shared_hom, round(sim, 4))


def run_task(task):
    """全局函数封装，用于多进程执行"""
    vcf_path, child, parent, mode = task
    child, parent, n_sites, n_shared, sim = calc_similarity(vcf_path, child, parent)
    return {
        "child_id": child,
        "parent_id": parent,
        "mode": mode,
        "n_sites": n_sites,
        "n_shared_hom": n_shared,
        "similarity": sim,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Calculate parent-child genotype similarity using homozygous sites"
    )
    parser.add_argument("--vcf", required=True, help="VCF file containing all samples")
    parser.add_argument("--children", required=True, help="Child sample list")
    parser.add_argument("--p1", help="Parent 1 sample list")
    parser.add_argument("--p2", help="Parent 2 sample list")
    parser.add_argument("--p_all", help="All parent samples (if not separated)")
    parser.add_argument("--out", required=True, help="Output TSV file")
    parser.add_argument("--threads", type=int, default=4, help="Parallel threads")
    args = parser.parse_args()

    children = parse_list(args.children)
    parents_1 = parse_list(args.p1)
    parents_2 = parse_list(args.p2)
    parents_all = parse_list(args.p_all)

    tasks = []
    if args.p_all:
        for c in children:
            for p in parents_all:
                tasks.append((args.vcf, c, p, "p_all"))
    else:
        for c in children:
            for p in parents_1:
                tasks.append((args.vcf, c, p, "p1"))
            for p in parents_2:
                tasks.append((args.vcf, c, p, "p2"))

    print(f"🧮 Total comparisons: {len(tasks)}")

    results = []
    with ProcessPoolExecutor(max_workers=args.threads) as ex:
        for res in ex.map(run_task, tasks):
            results.append(res)

    df = pd.DataFrame(results)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"✅ Result saved to {args.out}")


if __name__ == "__main__":
    main()
