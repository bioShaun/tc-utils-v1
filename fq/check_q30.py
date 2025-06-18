#!/usr/bin/env python3
"""
fastp Q30异常检测脚本
检测fastp输出中质量分数异常下降的位置
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np


def load_fastp_json(json_file):
    """加载fastp的JSON报告文件"""
    try:
        with open(json_file, "r") as f:
            data = json.load(f)
        return data
    except FileNotFoundError:
        print(f"错误: 找不到文件 {json_file}")
        sys.exit(1)
    except json.JSONDecodeError:
        print(f"错误: {json_file} 不是有效的JSON文件")
        sys.exit(1)


def extract_quality_data(fastp_data):
    """从fastp数据中提取质量信息"""
    quality_data = {}

    # 提取read1和read2的质量数据
    for read_type in ["read1", "read2"]:
        read_type_before_filter = f"{read_type}_before_filtering"
        if read_type in fastp_data:
            read_data = fastp_data[read_type_before_filter]
            if "quality_curves" in read_data:
                quality_data[read_type] = read_data["quality_curves"]

    return quality_data


def detect_q30_anomalies(
    read_type: str,
    quality_curves: Dict[str, List[float]],
    threshold=30,
    drop_threshold=0.1,
):
    """
    检测Q30异常

    参数:
    - quality_curves: 质量曲线数据
    - threshold: Q30阈值 (默认30)
    - drop_threshold: 质量下降阈值，基于头尾差异比例 (默认0.1，即10%)
    """
    anomalies = []

    for base in ["mean"]:
        if base in quality_curves:
            qualities = quality_curves[base]
            print(qualities)
            positions = list(range(len(qualities)))

            if len(qualities) < 10:  # 序列太短，跳过
                continue

            # 计算头部和尾部的平均质量
            head_size = min(10, len(qualities) // 4)  # 头部10个位置或1/4长度
            tail_size = min(10, len(qualities) // 4)  # 尾部10个位置或1/4长度

            head_quality = np.mean(qualities[:head_size])
            tail_quality = np.mean(qualities[-tail_size:])

            # 计算头尾差异
            head_tail_diff = head_quality - tail_quality

            # 基于头尾差异计算动态阈值
            if head_tail_diff > 0:
                dynamic_threshold = head_tail_diff * drop_threshold
            else:
                dynamic_threshold = 2.0  # 如果头部质量不比尾部高，使用固定阈值

            # 检测异常下降
            for i in range(1, len(qualities)):
                current_q = qualities[i]
                prev_q = qualities[i - 1]

                # 检测显著下降
                quality_drop = abs(prev_q - current_q)
                print(i - 1, i, quality_drop)
                if quality_drop > dynamic_threshold:
                    anomalies.append(
                        {
                            "base": base,
                            "read_type": read_type,
                            "position": positions[i],
                            "quality_before": prev_q,
                            "quality_after": current_q,
                            "drop": quality_drop,
                            "below_q30": current_q < threshold,
                            "head_quality": head_quality,
                            "tail_quality": tail_quality,
                            "head_tail_diff": head_tail_diff,
                            "dynamic_threshold": dynamic_threshold,
                        }
                    )

    return anomalies


def calculate_q30_rate(quality_curves, threshold=30):
    """计算Q30比率"""
    q30_rates = {}

    # 使用mean质量曲线计算Q30比率
    if "mean" in quality_curves:
        qualities = quality_curves["mean"]
        q30_count = sum(1 for q in qualities if q >= threshold)
        q30_rate = q30_count / len(qualities) if qualities else 0
        q30_rates["mean"] = q30_rate

    # 同时也计算各个碱基的Q30比率（如果存在）
    for base in ["A", "T", "C", "G"]:
        if base in quality_curves:
            qualities = quality_curves[base]
            q30_count = sum(1 for q in qualities if q >= threshold)
            q30_rate = q30_count / len(qualities) if qualities else 0
            q30_rates[base] = q30_rate

    return q30_rates


def plot_quality_curves(quality_data, output_dir="."):
    """绘制质量曲线图"""
    for read_type, quality_curves in quality_data.items():
        plt.figure(figsize=(12, 8))

        colors = {"A": "brown", "T": "red", "C": "blue", "G": "green", "mean": "black"}

        # 优先绘制mean曲线
        if "mean" in quality_curves:
            qualities = quality_curves["mean"]
            positions = list(range(len(qualities)))
            plt.plot(
                positions, qualities, color=colors["mean"], label="mean", linewidth=3
            )

        # 绘制各个碱基的曲线
        for base in ["A", "T", "C", "G"]:
            if base in quality_curves:
                qualities = quality_curves[base]
                positions = list(range(len(qualities)))
                plt.plot(
                    positions,
                    qualities,
                    color=colors[base],
                    label=base,
                    linewidth=2,
                    alpha=0.7,
                )

        # 添加Q30参考线
        plt.axhline(y=30, color="red", linestyle="--", alpha=0.7, label="Q30")

        plt.xlabel("Position")
        plt.ylabel("Quality Score")
        plt.title(f"Quality Curves - {read_type.upper()}")
        plt.legend()
        plt.grid(True, alpha=0.3)

        # 保存图片
        output_file = Path(output_dir) / f"quality_curves_{read_type}.png"
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"质量曲线图已保存: {output_file}")


def generate_report(fastp_data, quality_data, anomalies, output_file):
    """生成检测报告"""
    with open(output_file, "w") as f:
        f.write("fastp Q30异常检测报告\n")
        f.write("=" * 50 + "\n\n")

        # 基本信息
        f.write("基本信息:\n")
        f.write("-" * 20 + "\n")
        if "summary" in fastp_data:
            summary = fastp_data["summary"]
            if "before_filtering" in summary:
                f.write(
                    f"过滤前总reads数: {summary['before_filtering']['total_reads']:,}\n"
                )
                f.write(
                    f"过滤前总碱基数: {summary['before_filtering']['total_bases']:,}\n"
                )
                f.write(
                    f"过滤前Q30比率: {summary['before_filtering']['q30_rate']:.2%}\n"
                )
            if "after_filtering" in summary:
                f.write(
                    f"过滤后总reads数: {summary['after_filtering']['total_reads']:,}\n"
                )
                f.write(
                    f"过滤后总碱基数: {summary['after_filtering']['total_bases']:,}\n"
                )
                f.write(
                    f"过滤后Q30比率: {summary['after_filtering']['q30_rate']:.2%}\n"
                )
        f.write("\n")

        # Q30比率统计
        f.write("各碱基Q30比率:\n")
        f.write("-" * 20 + "\n")
        for read_type, quality_curves in quality_data.items():
            f.write(f"{read_type.upper()}:\n")
            q30_rates = calculate_q30_rate(quality_curves)
            for base, rate in q30_rates.items():
                f.write(f"  {base}: {rate:.2%}\n")
        f.write("\n")

        # 异常检测结果
        f.write("异常检测结果:\n")
        f.write("-" * 20 + "\n")
        if anomalies:
            f.write(f"发现 {len(anomalies)} 个异常位置:\n\n")
            for i, anomaly in enumerate(anomalies, 1):
                if i == 0:
                    f.write(
                        f"   头部平均质量: {anomaly['head_quality']:.2f}, "
                        f"尾部平均质量: {anomaly['tail_quality']:.2f}\n"
                    )
                    f.write(
                        f"   头尾差异: {anomaly['head_tail_diff']:.2f}, "
                        f"动态阈值: {anomaly['dynamic_threshold']:.2f}\n\n"
                    )
                    f.write("-" * 20 + "\n")
                f.write(
                    f"{anomaly['read_type']}: "
                    f"{i}. 位置: {anomaly['position']}, "
                    f"质量下降: {anomaly['drop']:.2f}, "
                    f"低于Q30: {'是' if anomaly['below_q30'] else '否'}\n"
                )
                f.write(
                    f"   质量变化: {anomaly['quality_before']:.2f} -> {anomaly['quality_after']:.2f}\n"
                )

        else:
            f.write("未发现显著的质量异常\n")

        # 建议
        f.write("建议:\n")
        f.write("-" * 20 + "\n")
        if anomalies:
            severe_anomalies = [a for a in anomalies if a["below_q30"]]
            if severe_anomalies:
                f.write("1. 发现质量严重下降的位置，建议:\n")
                f.write("   - 检查测序仪状态\n")
                f.write("   - 考虑增加质量过滤参数\n")
                f.write("   - 检查样品质量\n")
            else:
                f.write("1. 质量下降不严重，可能的原因:\n")
                f.write("   - 正常的末端质量下降\n")
                f.write("   - 轻微的测序质量波动\n")
        else:
            f.write("1. 质量控制良好，无需额外处理\n")


def main():
    parser = argparse.ArgumentParser(description="检测fastp输出中的Q30异常")
    parser.add_argument("json_file", help="fastp输出的JSON文件路径")
    parser.add_argument(
        "-o",
        "--output",
        default="q30_anomaly_report.txt",
        help="输出报告文件名 (默认: q30_anomaly_report.txt)",
    )
    parser.add_argument(
        "-t", "--threshold", type=float, default=30, help="Q30阈值 (默认: 30)"
    )
    parser.add_argument(
        "-d",
        "--drop-threshold",
        type=float,
        default=0.1,
        help="质量下降检测阈值，基于头尾差异的比例 (默认: 0.1，即10%%)",
    )
    parser.add_argument("--plot", action="store_true", help="生成质量曲线图")
    parser.add_argument("--plot-dir", default=".", help="图片输出目录 (默认: 当前目录)")

    args = parser.parse_args()

    # 加载fastp数据
    print(f"加载fastp数据: {args.json_file}")
    fastp_data = load_fastp_json(args.json_file)

    # 提取质量数据
    quality_data = extract_quality_data(fastp_data)
    if not quality_data:
        print("错误: 无法从JSON文件中提取质量数据")
        sys.exit(1)

    # 检测异常
    print("检测Q30异常...")
    all_anomalies = []
    for read_type, quality_curves in quality_data.items():
        anomalies = detect_q30_anomalies(
            read_type,
            quality_curves,
            threshold=args.threshold,
            drop_threshold=args.drop_threshold,
        )
        for anomaly in anomalies:
            anomaly["read_type"] = read_type
        all_anomalies.extend(anomalies)

    # 生成报告
    print(f"生成报告: {args.output}")
    generate_report(fastp_data, quality_data, all_anomalies, args.output)

    # 生成图片
    if args.plot:
        print("生成质量曲线图...")
        plot_quality_curves(quality_data, args.plot_dir)

    # 输出简要结果
    print(f"\n检测完成!")
    print(f"发现 {len(all_anomalies)} 个异常位置")
    if all_anomalies:
        severe_count = sum(1 for a in all_anomalies if a["below_q30"])
        print(f"其中 {severe_count} 个位置质量低于Q30")
    print(f"详细报告已保存至: {args.output}")


if __name__ == "__main__":
    main()
