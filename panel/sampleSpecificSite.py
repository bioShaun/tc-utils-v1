from pathlib import Path
import typer

import pandas as pd
from tqdm import tqdm


app = typer.Typer()


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


@app.command()
def gt_stats(
    gt_file: Path,
    all_sample_path: Path,
    case_sample_path: Path,
    case_name: str,
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
    gt_stats_file = gt_file.with_suffix(".stats")
    for i, gt_df in tqdm(enumerate(gt_df_list)):
        case_df = gt_df[case_sample_list]
        case_stats_df = get_gt_stats(case_df, case_name)
        control_name = f"non_{case_name}"
        control_df = gt_df[control_sample_list]
        control_stats_df = get_gt_stats(control_df, control_name)
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
    filter_columns = ["alt"]
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
