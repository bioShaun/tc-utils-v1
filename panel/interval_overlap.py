from pathlib import Path

import typer


def process_intervals(lines):
    # 按染色体号分组
    chrom_intervals = {}
    for line in lines:
        parts = line.strip().split()
        chrom = parts[0]
        label = f"{chrom}:{parts[1]}-{parts[2]}"
        if chrom not in chrom_intervals:
            chrom_intervals[chrom] = []
        chrom_intervals[chrom].append((int(parts[1]), int(parts[2]), label))

    result = []
    # 对每个染色体单独处理
    for chrom, intervals in chrom_intervals.items():
        # 收集当前染色体的所有分割点
        points = set()
        for start, end, _ in intervals:
            points.add(start)
            points.add(end)
        points = sorted(list(points))

        # 处理每个子区间
        for i in range(len(points) - 1):
            start = points[i]
            end = points[i + 1]
            # 找出覆盖当前区间的所有标签
            labels = []
            for interval_start, interval_end, label in intervals:
                if interval_start <= start and interval_end >= end:
                    labels.append(label)
            # 只有当有标签时才添加区间
            if labels:
                result.append((chrom, start, end, ",".join(sorted(labels))))

    # 按染色体号和起始位置排序
    result.sort(key=lambda x: (x[0], x[1]))
    return result


def main(bed: Path, out_file: Path):
    with open(bed, "r", encoding="utf-8") as f:
        lines = f.readlines()

    result = process_intervals(lines)
    with open(out_file, "w", encoding="utf-8") as f:
        for item in result:
            f.write("\t".join(map(str, item)) + "\n")


if __name__ == "__main__":
    typer.run(main)
