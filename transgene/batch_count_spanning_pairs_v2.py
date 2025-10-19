#!/usr/bin/env python3
"""
批量统计多个 BAM 文件中跨越插入边界 (Spanning Reads) 与插入内部 (Internal Reads) 的配对 reads，
计算 Spanning Ratio = Nspanning / (Nspanning + 0.5 * Ninternal)，
并自动判定样本基因型 (WT / Hetero / Homo)，输出精简表格。
"""

import argparse
from pathlib import Path

import pysam


def classify_genotype(spanning_ratio: float) -> str:
    """根据 Spanning Ratio 自动判定样本基因型"""
    if spanning_ratio >= 0.8:
        return "WT"  # 野生型，无插入
    elif spanning_ratio <= 0.2:
        return "HOM"  # 纯合插入
    else:
        return "HET"  # 杂合插入


def process_bam(bam_path, chrom, start, end, flank=10, mapq=0, outdir="."):
    """统计单个 BAM 文件中的 spanning/internal reads 并导出 SAM"""
    S, E, F = start, end, flank
    bam = pysam.AlignmentFile(bam_path, "rb")

    spanning_ids, internal_ids, total_ids = set(), set(), set()
    checked_reads = set()

    for region in [(S - F, S + F), (E - F, E + F)]:
        for aln in bam.fetch(chrom, region[0], region[1]):
            qn = aln.query_name
            if qn in checked_reads:
                continue
            checked_reads.add(qn)

            if aln.is_unmapped or aln.mate_is_unmapped:
                continue
            if aln.mapping_quality < mapq:
                continue
            if aln.next_reference_name not in (chrom, "="):
                continue

            rpos = aln.reference_start + 1
            rend = aln.reference_end
            mpos = aln.next_reference_start + 1
            total_ids.add(qn)

            # 一端在左，一端在右 → spanning
            if (rpos < S and mpos > E) or (rpos > E and mpos < S):
                spanning_ids.add(qn)
            # 两端都在插入区内 → internal
            elif (S <= rpos <= E) and (S <= mpos <= E):
                internal_ids.add(qn)

    bam.close()

    total = len(total_ids)
    n_spanning = len(spanning_ids)
    n_internal = len(internal_ids)
    denom = n_spanning + 0.5 * n_internal
    ratio_spanning = n_spanning / denom if denom > 0 else 0.0
    genotype = classify_genotype(ratio_spanning)

    # 输出 SAM 文件到子目录
    sample = Path(bam_path).stem
    sam_dir = Path(outdir) / "sam"
    sam_dir.mkdir(parents=True, exist_ok=True)

    out_spanning = sam_dir / f"{sample}_spanning.sam"
    out_internal = sam_dir / f"{sample}_internal.sam"

    bam = pysam.AlignmentFile(bam_path, "rb")
    span_sam = pysam.AlignmentFile(out_spanning, "w", header=bam.header)
    int_sam = pysam.AlignmentFile(out_internal, "w", header=bam.header)

    for aln in bam.fetch(until_eof=True):
        if aln.query_name in spanning_ids:
            span_sam.write(aln)
        elif aln.query_name in internal_ids:
            int_sam.write(aln)

    bam.close()
    span_sam.close()
    int_sam.close()

    return sample, n_spanning, n_internal, ratio_spanning, genotype


def main():
    ap = argparse.ArgumentParser(
        description="统计转基因插入区的 Spanning/Internal Reads 并判定基因型"
    )
    ap.add_argument(
        "--bams", nargs="+", required=True, help="输入 BAM 文件（支持通配符）"
    )
    ap.add_argument("--chrom", required=True, help="染色体名")
    ap.add_argument("--start", type=int, required=True, help="插入区起点")
    ap.add_argument("--end", type=int, required=True, help="插入区终点")
    ap.add_argument("--flank", type=int, default=10, help="边界检测范围 ±bp (默认10)")
    ap.add_argument("--mapq", type=int, default=0, help="最小 MAPQ (默认0)")
    ap.add_argument("--outdir", default="insertion_analysis", help="输出目录")
    args = ap.parse_args()

    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    (Path(args.outdir) / "sam").mkdir(parents=True, exist_ok=True)

    summary_path = Path(args.outdir) / "Transgene_Insertion_Summary.tsv"

    with open(summary_path, "w") as summary:
        # 精简后的列
        summary.write(
            "Sample\tChromosome\tInsertion_Start\tInsertion_End\t"
            "Spanning_Reads\tInternal_Reads\tSpanning_Ratio\tGenotype_Call\n"
        )

        for bam in args.bams:
            sample, n_span, n_int, ratio, genotype = process_bam(
                bam,
                args.chrom,
                args.start,
                args.end,
                args.flank,
                args.mapq,
                args.outdir,
            )

            summary.write(
                f"{sample}\t{args.chrom}\t{args.start}\t{args.end}\t"
                f"{n_span}\t{n_int}\t{ratio:.4f}\t{genotype}\n"
            )

            print(
                f"✅ {sample}: spanning={n_span}, internal={n_int}, ratio={ratio:.3f}, genotype={genotype}"
            )

    print(f"\n📄 Summary table written to: {summary_path}")
    print(f"📂 SAM files stored in: {Path(args.outdir) / 'sam'}")


if __name__ == "__main__":
    main()
