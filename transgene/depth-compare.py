import pandas as pd
import typer
from pathlib import Path
from typing import Annotated, List, Dict, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm  # 需要安装: pip install tqdm


def read_depth_median(depth_file: Path) -> Optional[float]:
    """
    读取depth文件并返回中位数。
    增加异常捕获，避免单个文件损坏导致程序中断。
    """
    if not depth_file.exists():
        return None

    try:
        # 显式指定 dtype 节省内存，bamdst depth 通常是整数，但在计算中位数时 pandas 会处理
        # 假设 header 存在且 Raw Depth 是第3列 (index 2)
        df = pd.read_csv(
            depth_file,
            sep="\t",
            compression="gzip",
            usecols=[2],
            dtype={2: "int32"},  # 深度通常为整数，用 int32 足够且省内存
        )
        # 检查是否读取到了数据
        if df.empty:
            return 0.0

        # 注意：如果文件有header，usecols=[2]会读取第三列。
        # 如果列名不是 "Raw Depth"，df.iloc[:, 0] 更通用
        return float(df.iloc[:, 0].median())

    except Exception as e:
        print(f"\nError reading {depth_file}: {e}")
        return None


def process_single_sample(
    transgene_dir: Path, background_bamdst_dir: Path
) -> Optional[Dict]:
    """
    处理单个样本的逻辑，方便放入并行池。
    """
    if not transgene_dir.is_dir():
        return None

    sample_id = transgene_dir.name
    trans_depth_file = transgene_dir / "depth.tsv.gz"
    bg_depth_file = background_bamdst_dir / sample_id / "depth.tsv.gz"

    # 检查背景文件是否存在，若不存在则跳过或返回特定标记
    if not bg_depth_file.exists():
        # 这里可以选择记录日志，或者返回 None
        return {
            "sample_id": sample_id,
            "transgene_depth": read_depth_median(trans_depth_file),
            "background_depth": None,
            "ratio": None,
            "note": "Background missing",
        }

    transgene_depth = read_depth_median(trans_depth_file)
    background_depth = read_depth_median(bg_depth_file)

    if transgene_depth is None or background_depth is None:
        return {
            "sample_id": sample_id,
            "transgene_depth": transgene_depth,
            "background_depth": background_depth,
            "ratio": 0.0,
            "note": "Read error",
        }

    ratio = transgene_depth / background_depth if background_depth > 0 else 0.0

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
    """比较转基因和背景BAM文件的深度（并行版）。"""

    # 获取所有待处理目录
    sample_dirs = [d for d in trans_gene_bamdst_dir.iterdir() if d.is_dir()]

    results = []

    print(f"Start processing {len(sample_dirs)} samples with {threads} threads...")

    # 使用多进程处理 (ProcessPoolExecutor 适合 CPU 密集型或 Pandas 操作)
    with ProcessPoolExecutor(max_workers=threads) as executor:
        # 提交任务
        futures = {
            executor.submit(process_single_sample, d, background_bamdst_dir): d.name
            for d in sample_dirs
        }

        # 使用 tqdm 显示进度条
        for future in tqdm(
            as_completed(futures), total=len(sample_dirs), unit="sample"
        ):
            res = future.result()
            if res:
                results.append(res)

    # 保存结果
    df_result = pd.DataFrame(results)

    # 调整列顺序，把 note 放在最后
    cols = ["sample_id", "transgene_depth", "background_depth", "ratio", "note"]
    # 确保列存在（防止全空情况）
    final_cols = [c for c in cols if c in df_result.columns]

    df_result = df_result[final_cols]
    df_result.to_csv(output_file, sep="\t", index=False)
    print(f"Results saved to {output_file}")


if __name__ == "__main__":
    typer.run(main)
