from pathlib import Path
from typing import List

import pandas as pd
import typer


def read_table_or_excel(file: Path) -> pd.DataFrame:
    """Read a file into a pandas DataFrame.

    The file can be a TSV or an Excel file.

    Args:
        file (Path): The file to read.

    Returns:
        pd.DataFrame: The DataFrame read from the file.

    Raises:
        ValueError: If the file's suffix is not .tsv or .xlsx.
    """
    read_functions = {".tsv": pd.read_table, ".xlsx": pd.read_excel}

    suffix = file.suffix
    if suffix in read_functions:
        return read_functions[suffix](file)

    raise ValueError(f"Only .tsv and .xlsx files are supported, not {suffix}")


def main(files: List[Path], output: Path) -> None:
    """Concatenate multiple TSV files into one.

    Args:
        files (List[Path]): List of files to concatenate.
        output (Path): Output file.
    """
    df = pd.concat([pd.read_table(file) for file in files], ignore_index=True)
    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)
