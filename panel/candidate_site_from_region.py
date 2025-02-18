from dataclasses import dataclass, field
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

    columns: List[str] = field(
        default_factory=lambda: ["sequence_score", "region_score", "maf"]
    )
    ascending: List[bool] = field(default_factory=lambda: [False, False, False])


def candidate_site_from_region(
    region_file: Path = typer.Argument(..., help="Path to the region file"),
    candidate_site_file: Path = typer.Argument(
        ..., help="Path to the candidate site file"
    ),
    output_file: Path = typer.Argument(..., help="Path to save the output"),
    max_region: int = typer.Option(
        100, "--max-region", "-m", help="Maximum number of regions to consider"
    ),
    candidate_per_region: int = typer.Option(
        1,
        "--candidate-per-region",
        "-c",
        help="Number of candidates to select per region",
    ),
    default_region_size: int = typer.Option(
        200, "--region-size", "-s", help="Default size of each region"
    ),
) -> None:
    """
    Select candidate sites from regions based on specified criteria.

    This tool processes region and candidate site files to select optimal candidate sites
    based on various parameters and scoring criteria.
    """

    region_df = pd.read_csv(
        region_file,
        sep="\t",
        header=None,
        names=["chr", "start", "end"],
        usecols=[0, 1, 2],
    )
    sites_df = pd.read_csv(candidate_site_file, sep="\t")
    sites_df.sort_values(
        by=PriorityOrder().columns, ascending=PriorityOrder().ascending, inplace=True
    )

    select_dfs = []
    for region in region_df.itertuples():
        start, end = int(region.start), int(region.end)
        span = end - start
        region_size = (
            default_region_size
            if default_region_size * max_region >= span
            else span // max_region
        )

        region_sites = sites_df[
            (sites_df["chrom"] == region.chr)
            & (sites_df["pos"] > start)
            & (sites_df["pos"] <= end)
        ].copy()

        region_sites["pos_range"] = pd.cut(
            region_sites["pos"], bins=range(start, end + 1, region_size), right=False
        )
        select_dfs.append(region_sites.groupby("pos_range").head(candidate_per_region))

    result_df = pd.concat(select_dfs)
    result_df.to_csv(output_file, sep="\t", index=False)
    typer.echo(f"Selected {len(result_df)} candidate sites.")


if __name__ == "__main__":
    typer.run(candidate_site_from_region)
