from dataclasses import dataclass
from pathlib import Path
from typing import List

import pandas as pd
import typer


@dataclass
class PriorityOrder:
    """
    PriorityOrder类用于定义优先级排序的规则。

    该类包含两个类属性：`columns`和`ascending`，分别指定了排序的列名和对应的排序顺序（升序或降序）。
    - `columns`: 包含排序依据的列名列表，默认为["sequence_score", "region_score", "maf"]。
    - `ascending`: 包含与`columns`对应的排序顺序列表，默认为[False, False, False]，表示所有列均按降序排序。

    该类主要用于数据排序，可以根据指定的列和顺序对数据进行优先级排序。

    示例：

    注意：
    - `columns`和`ascending`的长度必须一致，否则排序行为未定义。
    - 该类不包含构造函数，因此无需提供构造函数参数说明。
    """

    columns: List[str] = ["sequence_score", "region_score", "maf"]
    ascending: List[bool] = [False, False, False]


def candidate_site_from_region(
    region_file: Path,
    candidate_site_file: Path,
    output_file: Path,
    max_region: int = 100,
    candidate_per_region: int = 1,
    default_region_size: int = 200,
) -> None:
    """
    从给定的区域文件和候选位点文件中筛选出符合条件的候选位点，并将结果保存到指定的输出文件中。

    :param region_file: 区域文件的路径。
    :param candidate_site_file: 候选位点文件的路径。
    :param output_file: 输出文件的路径。
    :param max_region: 每个区域的最大分段数，默认值为100。
    :param candidate_per_region: 每个分段中选择的候选位点数，默认值为1。
    :param default_region_size: 默认的区域大小，默认值为200。
    """
    # 读取区域文件和候选位点文件
    region_df = pd.read_csv(
        region_file,
        sep="\t",
        header=None,
        names=["chr", "start", "end"],
        usecols=[0, 1, 2],
    )
    candidate_site_df = pd.read_csv(candidate_site_file, sep="\t")

    # 对候选位点进行排序
    candidate_site_df.sort_values(
        by=PriorityOrder.columns, ascending=PriorityOrder.ascending, inplace=True
    )

    # 初始化空列表用于存储筛选后的候选位点
    select_dfs = []

    # 遍历每个区域
    for region in region_df.itertuples():
        start = int(region.start)
        end = int(region.end)
        span = end - start
        region_size = (
            default_region_size
            if default_region_size * max_region >= span
            else span // max_region
        )

        # 选择当前区域内的位点
        region_df = candidate_site_df[
            (candidate_site_df["chrom"] == region.chrom)
            & (candidate_site_df["pos"] > start)
            & (candidate_site_df["pos"] <= end)
        ].copy()

        # 按分段分组并选择位点
        region_df["pos_range"] = pd.cut(
            region_df["pos"], bins=range(start, end + 1, region_size), right=False
        )
        region_best_df = region_df.groupby("pos_range").head(candidate_per_region)

        # 添加到列表中
        select_dfs.append(region_best_df)

    # 合并筛选后的位点并保存到输出文件
    select_df = pd.concat(select_dfs)
    select_df.to_csv(output_file, sep="\t", index=False)

    # 输出筛选后的候选位点数量
    typer.echo(f"Selected {len(select_df)} candidate sites.")


if __name__ == "__main__":
    typer.run(candidate_site_from_region)
