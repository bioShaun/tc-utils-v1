#!/usr/bin/env python3
"""
fastp Q30异常检测脚本
检测fastp输出中质量分数异常下降的位置
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm


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
        if read_type_before_filter in fastp_data:
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


def generate_r1_r2_info(read_type: Literal["read1", "read2"], anomalies: List[Dict]):
    data_info = [each for each in anomalies if each["read_type"] == read_type]
    if data_info:
        bad_number = len(data_info)
        bad_start = data_info[0]["position"]
        bad_end = data_info[-1]["position"]
        max_drop = max(each["drop"] for each in data_info)
    else:
        bad_number = 0
        bad_start = 0
        bad_end = 0
        max_drop = 0
    return {
        f"{read_type}异常位点数量": bad_number,
        f"{read_type}异常位点起始": bad_start,
        f"{read_type}异常位点结束": bad_end,
        f"{read_type}最大质量下降": max_drop,
    }


def generate_report_data(anomalies, sample_name):
    """生成检测报告"""
    r1_info = generate_r1_r2_info("read1", anomalies)
    r2_info = generate_r1_r2_info("read2", anomalies)
    return {
        "libid": sample_name,
        **r1_info,
        **r2_info,
        "Total异常位点数量": r1_info["read1异常位点数量"]
        + r2_info["read2异常位点数量"],
    }


def main():
    parser = argparse.ArgumentParser(description="检测fastp输出中的Q30异常")
    parser.add_argument("json_dir", help="fastp输出的JSON文件路径目录")
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
    json_file_path_list = list(Path(args.json_dir).glob("*.json"))
    report_data_list = []
    for json_file_path in tqdm(json_file_path_list):
        print(f"加载fastp数据: {json_file_path}")
        fastp_data = load_fastp_json(json_file_path)

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
        sample_name = Path(json_file_path).stem
        report_data = generate_report_data(all_anomalies, sample_name)
        report_data_list.append(report_data)
    report_data_df = pd.DataFrame(report_data_list)
    report_data_df.to_csv(args.output, sep="\t", index=False, float_format="%.3f")


if __name__ == "__main__":
    main()
