from pathlib import Path
from typing import Tuple

import pandas as pd
import typer


def extract_seq(seq: str) -> Tuple[int, str]:
    left = seq.split("[")[0]
    right = seq.split("]")[1]
    center = seq.split("/")[0][-1]
    return len(left), left + center + right


def main(axiom_annotation: Path) -> None:
    if axiom_annotation.suffix == ".csv":
        axiom_df = pd.read_csv(axiom_annotation, comment="#")
    elif axiom_annotation.suffix == ".xlsx":
        axiom_df = pd.read_excel(axiom_annotation, comment="#")
    else:
        raise ValueError("Unknown file format")
    seq_info = axiom_df["Flank"].apply(extract_seq)
    axiom_df["sequence"] = [each[1] for each in seq_info]
    axiom_df["allele_position"] = [each[0] for each in seq_info]
    for _, row in axiom_df.iterrows():
        print(f">{row['Probe Set ID']}_{row['allele_position']}\n{row['sequence']}")


if __name__ == "__main__":
    typer.run(main)
