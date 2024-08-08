from pathlib import Path

import pandas as pd
import typer


def tsv2excel(tsv_path: Path, excel_path: Path):
    df = pd.read_table(tsv_path)
    df.to_excel(excel_path, index=False)


if __name__ == "__main__":
    typer.run(tsv2excel)
