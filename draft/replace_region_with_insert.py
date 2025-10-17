#!/usr/bin/env python3
"""
使用 pyfaidx 将转基因序列插入到基因组指定区域，
并输出完整基因组（所有染色体）。
"""

import argparse

from pyfaidx import Fasta


def load_insert_sequence(insert_fasta: str) -> str:
    """读取插入片段序列（假设只有一个序列）"""
    insert = Fasta(insert_fasta)
    if len(insert.keys()) != 1:
        raise ValueError("插入序列 FASTA 应只包含一个序列")
    name = list(insert.keys())[0]
    return str(insert[name][:].seq)


def replace_sequence(
    ref_fa: str, chrom: str, start: int, end: int, insert_seq: str, out_fa: str
):
    """
    将 ref_fa 中指定区域 [start, end] 替换为 insert_seq，
    输出完整基因组。
    注意: 坐标为 1-based 闭区间。
    """
    ref = Fasta(ref_fa)
    if chrom not in ref.keys():
        raise ValueError(f"参考基因组中不存在染色体 {chrom}")

    with open(out_fa, "w") as out:
        for chr_name in ref.keys():
            seq = ref[chr_name][:].seq

            if chr_name == chrom:
                print(f"🧬 替换 {chrom}:{start}-{end} 区域...")
                left = seq[: start - 1]
                right = seq[end:]
                seq = left + insert_seq + right
                print(f"✅ 已插入 {len(insert_seq)} bp 新序列")

            # 按 60 bp 一行写入
            out.write(f">{chr_name}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i : i + 60] + "\n")

    print(f"\n🎉 替换完成，输出文件: {out_fa}")


def main():
    parser = argparse.ArgumentParser(
        description="将转基因序列插入到指定基因组区域（保留所有染色体）"
    )
    parser.add_argument("--ref", required=True, help="参考基因组 fasta 文件")
    parser.add_argument("--insert", required=True, help="插入序列 fasta 文件")
    parser.add_argument("--chrom", required=True, help="染色体名称")
    parser.add_argument(
        "--start", required=True, type=int, help="插入起始位置 (1-based)"
    )
    parser.add_argument("--end", required=True, type=int, help="插入结束位置 (1-based)")
    parser.add_argument("--out", required=True, help="输出 fasta 文件")
    args = parser.parse_args()

    insert_seq = load_insert_sequence(args.insert)
    replace_sequence(args.ref, args.chrom, args.start, args.end, insert_seq, args.out)


if __name__ == "__main__":
    main()
