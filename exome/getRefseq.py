import re
import typer

import pandas as pd


from pathlib import Path


def extract_ref(sequence: str) -> str:
    match_res = re.match("(\w+)\[(\w+)/\w+\](\w+)", sequence)
    if match_res:
        return "".join(match_res.groups())
    else:
        return ""


def main(snp_file: Path) -> None:
    df = pd.read_csv(snp_file)
    df["ref_seq"] = df["Sequence"].map(extract_ref)
    for _, row in df.iterrows():
        print(f">{row['Name']}")
        print(row["ref_seq"])


if __name__ == "__main__":
    typer.run(main)
