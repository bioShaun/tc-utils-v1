from pathlib import Path
from typing import Tuple

import pandas as pd
import typer
from tqdm import tqdm

app = typer.Typer()


GT_COLUMN_PREFIX = ["chrom", "pos", "ref", "alt", "filter"]


GT_MAP = {
    "./.": "miss",
    "0/0": "ref",
    "0/1": "het",
    "1/1": "alt",
}


def get_gt_stats(df: pd.DataFrame, name: str) -> pd.DataFrame:
    stats_df = pd.DataFrame(df.apply(lambda x: x.value_counts() / df.shape[1], axis=1))  # type: ignore
    gt_list = list(GT_MAP.keys())
    for gt in gt_list:
        if gt not in stats_df.columns:
            stats_df[gt] = 0
    stats_df = stats_df[gt_list]
    stats_df.columns = [f"{name}_{GT_MAP[each]}" for each in stats_df.columns]
    stats_df.fillna(0, inplace=True)
    return stats_df


def get_gt_type(gt: str) -> str:
    if gt in ["./.", "."]:
        return "miss"
    try:
        allele1, allele2 = gt.split("/")[:2]
    except ValueError as exc:
        print("error gt:", gt)
        raise ValueError(f"{gt} format error!") from exc
    if allele1 != allele2:
        return "het"
    if allele1 == "0":
        return "ref"
    return "alt"


def allele_stats(row: pd.Series) -> Tuple[float, float, float]:
    gt_type = row.map(get_gt_type)
    allele_count = gt_type.value_counts()
    miss_count = allele_count.get("miss", 0)
    het_count = allele_count.get("het", 0)
    real_sample_count = len(row) - miss_count
    miss_rate = miss_count / len(row)
    het_rate = het_count / real_sample_count
    alt_count = allele_count.get("alt", 0)
    alt_rate = (alt_count * 2 + het_count) / (real_sample_count * 2)
    return miss_rate, het_rate, min(alt_rate, 1 - alt_rate)


@app.command()
def gt_stats(
    gt_file: Path,
    all_sample_path: Path,
    case_sample_path: Path,
    control_sample_path: Path,
    case_name: str,
    control_name: str,
    test: bool = False,
    chunck_size: int = 10_000,
):
    all_sample_list = [each.strip() for each in all_sample_path.open()]
    case_sample_list = [each.strip() for each in case_sample_path.open()]
    control_sample_list = [each.strip() for each in control_sample_path.open()]
    if test:
        gt_df_list = pd.read_table(
            gt_file,
            header=None,
            names=[*GT_COLUMN_PREFIX, *all_sample_list],
            index_col=[0, 1, 2, 3, 4],
            nrows=100,
            chunksize=chunck_size,
        )
    else:
        gt_df_list = pd.read_table(
            gt_file,
            header=None,
            names=[*GT_COLUMN_PREFIX, *all_sample_list],
            index_col=[0, 1, 2, 3, 4],
            chunksize=chunck_size,
        )
    gt_stats_file = gt_file.with_suffix(".stats")
    for i, gt_df in tqdm(enumerate(gt_df_list)):
        case_df = gt_df[case_sample_list]
        case_stats_df = pd.DataFrame(
            case_df.apply(allele_stats, axis=1), columns=["stats_info"]
        )
        case_stats_df[
            [f"{case_name}_missing", f"{case_name}_het", f"{case_name}_maf"]
        ] = pd.DataFrame(
            case_stats_df["stats_info"].tolist(), index=case_stats_df.index
        )
        case_stats_df.drop("stats_info", axis=1, inplace=True)
        control_df = gt_df[control_sample_list]
        control_stats_df = pd.DataFrame(
            control_df.apply(allele_stats, axis=1), columns=["stats_info"]
        )
        control_stats_df[
            [f"{control_name}_missing", f"{control_name}_het", f"{control_name}_maf"]
        ] = pd.DataFrame(
            control_stats_df["stats_info"].tolist(), index=case_stats_df.index
        )
        control_stats_df.drop("stats_info", axis=1, inplace=True)
        merged_df = case_stats_df.merge(
            control_stats_df, left_index=True, right_index=True
        )
        if i == 0:
            header = True
            mode = "w"
        else:
            header = False
            mode = "a"
        merged_df.to_csv(gt_stats_file, sep="\t", mode=mode, header=header)


@app.command()
def gt_stats_filter(
    stats_file: Path,
    case_name: str,
    case_va_portion: float = 0.8,
    contral_va_portion: float = 0.01,
    allow_het: bool = True,
):
    stats_df = pd.read_table(stats_file)
    filter_columns = ["alt", "ref"]
    if allow_het:
        filter_columns.append("het")
    for col in filter_columns:
        case_col = f"{case_name}_{col}"
        control_col = f"non_{case_name}_{col}"
        filter1 = stats_df[case_col] >= case_va_portion
        filter2 = stats_df[control_col] < contral_va_portion
        filter_df = stats_df[filter1 & filter2]
        filter_file = stats_file.with_suffix(f".filter_{col}.tsv")
        filter_df.to_csv(filter_file, sep="\t", index=False)


if __name__ == "__main__":
    app()
