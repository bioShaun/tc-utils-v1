from enum import Enum
from pathlib import Path

import pandas as pd
import typer


class OutType(str, Enum):
    txt = "txt"
    xlsx = "xlsx"


LOCATION_COLS = ["CHROM", "POS", "REF", "ALT"]


def main(dp_file: Path, sample_file: Path, out_file: Path) -> None:
    sample_list = pd.read_csv(sample_file, header=None)[0].tolist()
    dp_df = pd.read_table(dp_file, header=None, names=[*LOCATION_COLS, *sample_list])
    dp_df = dp_df.set_index(LOCATION_COLS)
    dp_df = dp_df.apply(lambda x: x if str(x).isdigit() else 0).astype("int")
    dp_df.to_excel(out_file)


if __name__ == "__main__":
    typer.run(main)
