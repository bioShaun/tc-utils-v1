import polars as pl
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
        # 使用polars读取第3列 (Raw Depth)
        # polars 可以自动识别 gzip 压缩
        df = pl.read_csv(
            depth_file,
            separator="\t",
        )

        # 获取 Raw Depth 列
        depth_series = df["Raw Depth"]

        # 过滤掉深度 <= 0 的点，仅保留被测序覆盖到的区域
        valid_depths = depth_series.filter(depth_series > 0)

        # 检查过滤后是否为空
        if valid_depths.is_empty():
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

    transgene_depth = read_depth_median(trans_depth_file)

    # 背景文件不存在时返回None
    if not bg_depth_file.exists():
        background_depth = None
    else:
        background_depth = read_depth_median(bg_depth_file)

    # 检查读取是否成功
    if transgene_depth is None or background_depth is None:
        return None

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

    # 按sample_id排序
    results.sort(key=lambda x: x["sample_id"])

    df_result = pl.DataFrame(results)

    # 输出
    df_result.write_csv(output_file, separator="\t")
    print(f"Results saved to {output_file}")


if __name__ == "__main__":
    typer.run(main)
