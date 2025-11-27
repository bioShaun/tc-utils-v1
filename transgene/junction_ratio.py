#!/usr/bin/env python3
"""
从 BAM 文件中提取指定区域 (±flank) 内的 reads，
并根据是否跨越插入区域，将它们分别输出为两个 SAM 文件。
"""

import argparse

import pysam


def split_junction_reads(
    bam_path, chrom, start, end, flank=2, mapq_cutoff=10, out_prefix="junction_split"
):
    """
    提取目标区域 reads 并分成两个 SAM 文件:
    - <prefix>_junction.sam
    - <prefix>_nonjunction.sam
    Junction 定义: read.reference_start < start 且 read.reference_end > end
    """
    bam = pysam.AlignmentFile(bam_path, "rb")

    region_start = start - flank
    region_end = end + flank
    total_reads = 0
    junction_reads = 0
    non_junction_reads = 0

    out_junc = out_prefix + "_junction.sam"
    out_non = out_prefix + "_nonjunction.sam"

    with pysam.AlignmentFile(
        out_junc, "w", header=bam.header
    ) as junc_sam, pysam.AlignmentFile(out_non, "w", header=bam.header) as non_sam:
        for read in bam.fetch(chrom, region_start, region_end):
            if read.is_unmapped or read.mapping_quality < mapq_cutoff:
                continue
            total_reads += 1

            ref_start = read.reference_start
            ref_end = read.reference_end
            print(read.reference_start, read.reference_end)

            if ref_start < start and ref_end > end:
                junc_sam.write(read)
                junction_reads += 1
            else:
                non_sam.write(read)
                non_junction_reads += 1

    bam.close()

    ratio = junction_reads / total_reads if total_reads > 0 else 0

    print(f"📍 Region: {chrom}:{start}-{end} ±{flank}bp")
    print(f"Total reads: {total_reads}")
    print(f"Junction reads: {junction_reads}")
    print(f"Non-junction reads: {non_junction_reads}")
    print(f"Junction ratio: {ratio:.4f}")
    print(f"✅ Junction SAM: {out_junc}")
    print(f"✅ Non-junction SAM: {out_non}")

    return {
        "total_reads": total_reads,
        "junction_reads": junction_reads,
        "non_junction_reads": non_junction_reads,
        "ratio": ratio,
        "junction_sam": out_junc,
        "non_junction_sam": out_non,
    }


def main():
    parser = argparse.ArgumentParser(
        description="分离 junction 与非-junction reads 并输出原始 SAM"
    )
    parser.add_argument("--bam", required=True, help="输入 BAM 文件")
    parser.add_argument("--chrom", required=True, help="染色体名称")
    parser.add_argument("--start", required=True, type=int, help="插入起始位置")
    parser.add_argument("--end", required=True, type=int, help="插入结束位置")
    parser.add_argument("--flank", type=int, default=2, help="上下游扩展范围(bp)")
    parser.add_argument("--mapq", type=int, default=10, help="最小比对质量过滤")
    parser.add_argument("--out", default="junction_split", help="输出文件前缀")
    args = parser.parse_args()

    split_junction_reads(
        bam_path=args.bam,
        chrom=args.chrom,
        start=args.start,
        end=args.end,
        flank=args.flank,
        mapq_cutoff=args.mapq,
        out_prefix=args.out,
    )


if __name__ == "__main__":
    main()
