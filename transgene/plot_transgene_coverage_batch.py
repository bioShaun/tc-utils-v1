#!/usr/bin/env python3
"""
批量绘制多个样品的转基因插入区 reads 覆盖密度图：
- 输入为 depth 目录 (每个样本一个子目录，内含 depth.tsv.gz)
- y轴从0开始
- 覆盖用颜色填充表示
- 标注插入起点与终点位置
- 输出命名格式: <sample>_coverage-density.(png|pdf)
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_coverage_colormap(depth_file, chrom, start, end, sample, outdir="."):
    """绘制单样本覆盖密度图"""
    df = pd.read_csv(depth_file, sep="\t", header=0, compression="infer")

    # 自动识别深度列
    depth_col = "Raw Depth"
    if depth_col not in df.columns:
        raise ValueError("❌ 未检测到包含 'depth' 的列，请检查输入文件。")

    # 筛选目标区域 ±50bp
    df = df[df[df.columns[0]].astype(str) == str(chrom)]
    if df.empty:
        print(f"⚠️ {sample}: 区域 {chrom}:{start}-{end} 无覆盖数据，跳过。")
        return

    # 绘图
    plt.figure(figsize=(10, 5))
    x = df["Pos"]
    y = df[depth_col]

    plt.fill_between(x, y, color="royalblue", alpha=0.5, label="Read Coverage")
    plt.axvline(
        start, color="red", linestyle="--", lw=1.3, label=f"Insert Start ({start})"
    )
    plt.axvline(
        end, color="orange", linestyle="--", lw=1.3, label=f"Insert End ({end})"
    )

    plt.ylim(bottom=0)
    plt.xlabel("Genomic Position (bp)")
    plt.ylabel("Read Depth")
    plt.title(f"{sample} Coverage around Transgene Insertion: {chrom}:{start}-{end}")

    plt.legend(frameon=False, fontsize=9)
    plt.grid(alpha=0.3, linestyle="--", linewidth=0.5)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_png = outdir / f"{sample}_coverage-density.png"
    out_pdf = outdir / f"{sample}_coverage-density.pdf"

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()

    print(f"✅ {sample}: 图已生成 -> {out_png}")


def main():
    ap = argparse.ArgumentParser(description="批量绘制转基因插入区域 reads 覆盖密度图")
    ap.add_argument(
        "--depth_dir", required=True, help="输入深度目录（每个样本一个子目录）"
    )
    ap.add_argument("--chrom", required=True, help="染色体名")
    ap.add_argument("--start", type=int, required=True, help="插入区起点")
    ap.add_argument("--end", type=int, required=True, help="插入区终点")
    ap.add_argument("--outdir", default="coverage_plot", help="输出目录")
    args = ap.parse_args()

    depth_dir = Path(args.depth_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 扫描所有样品子目录
    depth_files = list(depth_dir.glob("*/depth.tsv*"))
    if not depth_files:
        raise FileNotFoundError(f"❌ 未在 {depth_dir} 下找到任何 depth.tsv 文件。")

    print(f"🔍 检测到 {len(depth_files)} 个样品，将绘制覆盖图...")
    for f in sorted(depth_files):
        sample = f.parent.name
        try:
            plot_coverage_colormap(f, args.chrom, args.start, args.end, sample, outdir)
        except Exception as e:
            print(f"❌ {sample}: 绘图失败 - {e}")

    print(f"\n📂 所有结果已输出到目录: {outdir.resolve()}")


if __name__ == "__main__":
    main()
