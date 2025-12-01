from pathlib import Path

import pandas as pd
import typer
from tqdm import tqdm

INFO_COLUMNS = [
    "Location",
    "Variant_type",
    "Chr",
    "Pos",
    "Ref",
    "Alt",
    "Gene_region",
    "Gene_symbol",
    "Anno",
]


def parse_sample_series(cols: list[str]) -> dict[str, list[str]]:
    group_dict = {}
    for col in cols:
        group_id = col.split("-")[0]
        group_dict.setdefault(group_id, []).append(col)
    return group_dict


def main(va_table: Path) -> None:
    va_df = pd.read_excel(va_table)
    sample_cols = [c for c in va_df.columns if c not in INFO_COLUMNS]
    sample_df = va_df[sample_cols]
    sample_group_dict = parse_sample_series(sample_cols)
    for group_id, group_cols in tqdm(sample_group_dict.items()):
        group_freq_cols = [c for c in group_cols if c.endswith("_Variant_Frequence")]
        group_freq_df = sample_df[group_freq_cols]
        filter_group_freq_df = group_freq_df > 0.001
        filter_group_freq_sum = filter_group_freq_df.sum(axis=1)
        passed_filter_group_freq_sum = filter_group_freq_sum[
            (filter_group_freq_sum >= 1) & (filter_group_freq_sum <= 3)
        ]
        group_out_cols = [*INFO_COLUMNS[:6], *group_cols, *INFO_COLUMNS[6:]]
        group_df = va_df.loc[passed_filter_group_freq_sum.index][group_out_cols]
        outdir = va_table.parent
        out_name = f"{group_id}.filter.{va_table.name}"
        group_filter_file = outdir / out_name
        group_df.to_excel(group_filter_file, index=False)


if __name__ == "__main__":
    typer.run(main)
