import typer
import pandas as pd

from pathlib import Path


def main(anno_table: Path, out_dir: Path, probe_name: Path) -> None:
    if anno_table.suffix == ".xlsx":
        anno_df = pd.read_excel(anno_table)
    else:
        anno_df = pd.read_table(anno_table)
    anno_df["pos_0"] = anno_df["pos"] - 1
