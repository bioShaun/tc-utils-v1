import pandas as pd
import typer
from pathlib import Path
from typing import Annotated, Dict, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm


def read_depth_median(depth_file: Path) -> Optional[float]:
    """
    读取depth文件，去除深度为0的点后返回中位数。
    """
    if not depth_file.exists():
        return None

    try:
        # 1. 读取数据
        # usecols=[2]: 只读取第3列 (Raw Depth)
        # dtype={2: "int32"}: 强制转换为int32以节省内存
        df = pd.read_csv(
            depth_file, sep="\t", compression="gzip", usecols=[2], dtype={2: "int32"}
        )

        if df.empty:
            return 0.0

        # 获取深度列数据 (Series)
        depth_series = df.iloc[:, 0]

        # 2. 关键修改：过滤掉深度 <= 0 的点
        # 仅保留被测序覆盖到的区域
        valid_depths = depth_series[depth_series > 0]

        # 3. 检查过滤后是否为空 (即原文件全是0的情况)
        if valid_depths.empty:
            return 0.0

        return float(valid_depths.median())

    except Exception as e:
        print(f"\nError reading {depth_file}: {e}")
        return None


def process_single_sample(
    transgene_dir: Path, background_bamdst_dir: Path
) -> Optional[Dict]:
    """处理单个样本的逻辑"""
    if not transgene_dir.is_dir():
        return None

    sample_id = transgene_dir.name
    trans_depth_file = transgene_dir / "depth.tsv.gz"
    bg_depth_file = background_bamdst_dir / sample_id / "depth.tsv.gz"

    # 检查背景文件
    if not bg_depth_file.exists():
        return {
            "sample_id": sample_id,
            "transgene_depth": read_depth_median(trans_depth_file),
            "background_depth": None,
            "ratio": None,
            "note": "Background missing",
        }

    # 读取两个文件的中位数
    transgene_depth = read_depth_median(trans_depth_file)
    background_depth = read_depth_median(bg_depth_file)

    # 检查读取是否成功
    if transgene_depth is None or background_depth is None:
        return {
            "sample_id": sample_id,
            "transgene_depth": transgene_depth,
            "background_depth": background_depth,
            "ratio": 0.0,
            "note": "Read error",
        }

    # 计算比率
    ratio = (
        transgene_depth / background_depth
        if background_depth and background_depth > 0
        else 0.0
    )

    return {
        "sample_id": sample_id,
        "transgene_depth": transgene_depth,
        "background_depth": background_depth,
        "ratio": ratio,
        "note": "OK",
    }


def main(
    trans_gene_bamdst_dir: Annotated[
        Path, typer.Argument(help="转基因 bamdst 结果目录")
    ],
    background_bamdst_dir: Annotated[Path, typer.Argument(help="背景 bamdst 结果目录")],
    output_file: Annotated[Path, typer.Argument(help="输出文件路径")],
    threads: Annotated[int, typer.Option(help="并行处理的线程数")] = 4,
) -> None:
    """
    比较转基因和背景BAM文件的深度。
    会自动去除 coverage 为 0 的位点计算中位数。
    """

    sample_dirs = [d for d in trans_gene_bamdst_dir.iterdir() if d.is_dir()]
    results = []

    print(
        f"Start processing {len(sample_dirs)} samples with {threads} threads (filtering 0-depth sites)..."
    )

    # 并行处理
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(process_single_sample, d, background_bamdst_dir): d.name
            for d in sample_dirs
        }

        for future in tqdm(
            as_completed(futures), total=len(sample_dirs), unit="sample"
        ):
            res = future.result()
            if res:
                results.append(res)

    # 保存结果
    if not results:
        print("No results generated.")
        return

    df_result = pd.DataFrame(results)

    # 整理列顺序
    cols = ["sample_id", "transgene_depth", "background_depth", "ratio", "note"]
    final_cols = [c for c in cols if c in df_result.columns]

    df_result = df_result[final_cols]
    df_result.to_csv(output_file, sep="\t", index=False)
    print(f"Results saved to {output_file}")


if __name__ == "__main__":
    typer.run(main)
