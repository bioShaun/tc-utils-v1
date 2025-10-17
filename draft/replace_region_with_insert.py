#!/usr/bin/env python3
"""
replace_region_with_insert.py
---------------------------------
使用 pyfaidx 将参考基因组中指定区间 (start–end) 替换为给定的插入序列。

示例：
    python replace_region_with_insert.py \
        --ref ref.fa \
        --insert insert.fa \
        --chrom chr6 \
        --start 150747133 \
        --end 150747159 \
        --out ref_with_insert.fa

参数说明：
    --ref     原始参考基因组FASTA
    --insert  要插入的序列FASTA
    --chrom   染色体名称（必须与ref中一致）
    --start   替换区间起点 (1-based)
    --end     替换区间终点 (1-based, 包含)
    --out     输出新的参考基因组FASTA
"""

import argparse
from pathlib import Path

from pyfaidx import Fasta


def read_insert_seq(insert_fa: str) -> str:
    """读取插入序列（取第一条FASTA记录）"""
    ins = Fasta(insert_fa)
    first = list(ins.keys())[0]
    seq = str(ins[first][:].seq).upper()
    ins.close()
    return seq


def replace_sequence(
    ref_fa: str, chrom: str, start: int, end: int, insert_seq: str, out_fa: str
):
    """将参考基因组指定区间替换为插入序列"""
    ref = Fasta(ref_fa, as_raw=True)
    if chrom not in ref:
        raise ValueError(f"染色体 {chrom} 不存在于参考基因组中")

    chrom_seq = ref[chrom][:].seq
    seq_len = len(chrom_seq)
    if start < 1 or end > seq_len or start > end:
        raise ValueError(f"无效区间: {start}-{end}, {chrom} 长度 {seq_len}")

    # 替换: start 和 end 是 1-based，且 end 包含
    left = chrom_seq[: start - 1]
    right = chrom_seq[end:]
    new_seq = left + insert_seq + right

    # 写出新的参考FASTA
    with open(out_fa, "w") as out:
        for name in ref.keys():
            if name == chrom:
                out.write(f">{chrom}\n")
                for i in range(0, len(new_seq), 60):
                    out.write(new_seq[i : i + 60] + "\n")
            else:
                seq = ref[name][:].seq
                out.write(f">{name}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i : i + 60] + "\n")

    ref.close()
    print(f"✅ 新参考文件已生成: {out_fa}")
    print(f"替换区域: {chrom}:{start}-{end} ({end - start + 1} bp)")
    print(f"插入序列长度: {len(insert_seq)} bp")
    print(f"新染色体长度变化: {len(insert_seq) - (end - start + 1)} bp")


def main():
    parser = argparse.ArgumentParser(description="将参考基因组指定区间替换为插入序列")
    parser.add_argument("--ref", required=True, help="参考基因组 FASTA 文件")
    parser.add_argument("--insert", required=True, help="插入序列 FASTA 文件")
    parser.add_argument("--chrom", required=True, help="染色体名称")
    parser.add_argument(
        "--start", required=True, type=int, help="起始位置 (1-based, 包含)"
    )
    parser.add_argument(
        "--end", required=True, type=int, help="终止位置 (1-based, 包含)"
    )
    parser.add_argument("--out", required=True, help="输出新参考基因组文件")
    args = parser.parse_args()

    insert_seq = read_insert_seq(args.insert)
    replace_sequence(args.ref, args.chrom, args.start, args.end, insert_seq, args.out)


if __name__ == "__main__":
    main()
