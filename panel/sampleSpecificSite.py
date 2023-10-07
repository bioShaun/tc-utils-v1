from pathlib import Path
import typer

import pandas as pd
from tqdm import tqdm


GT_COLUMN_PREFIX = ["chrom", "pos", "ref", "alt", "filter"]


GT_MAP = {
    "./.": "miss",
    "0/0": "ref",
    "1/1": "alt",
    "0/1": "het",
}


def get_gt_stats(df: pd.DataFrame, name: str) -> pd.DataFrame:
    stats_df = pd.DataFrame(df.apply(lambda x: x.value_counts() / df.shape[1], axis=1))  # type: ignore
    stats_df.columns = [f"{name}_{GT_MAP[each]}" for each in stats_df.columns]
    stats_df.fillna(0, inplace=True)
    return stats_df


def main(
    gt_file: Path,
    all_sample_path: Path,
    case_sample_path: Path,
    case_name: str,
    case_va_portion: float = 0.9,
    contral_va_portion: float = 0.05,
    control_va_count: int = 10,
    test: bool = False,
    chunck_size: int = 10_000,
):
    all_sample_list = [each.strip() for each in all_sample_path.open()]
    case_sample_list = [each.strip() for each in case_sample_path.open()]
    control_sample_list = [
        each for each in all_sample_list if each not in case_sample_list
    ]
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
    stats_df_list = []
    for gt_df in tqdm(gt_df_list):
        case_df = gt_df[case_sample_list]
        case_stats_df = get_gt_stats(case_df, case_name)
        control_name = f"non_{case_name}"
        control_df = gt_df[control_sample_list]
        control_stats_df = get_gt_stats(control_df, control_name)
        merged_df = case_stats_df.merge(
            control_stats_df, left_index=True, right_index=True
        )
        stats_df_list.append(merged_df)
    stats_df = pd.concat(stats_df_list)
    gt_stats_file = gt_file.with_suffix(".stats")
    stats_df.to_csv(gt_stats_file, sep="\t")


if __name__ == "__main__":
    typer.run(main)
