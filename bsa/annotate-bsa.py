import pandas as pd
import typer
from pathlib import Path

import shutil
import os
from tqdm import tqdm


def copy_directory(source_dir: Path, destination_dir: Path, overwrite=False) -> bool:
    """
    Copy a directory from source to destination.

    Args:
        source_dir (str): Source directory path
        destination_dir (str): Destination directory path
        overwrite (bool): If True, overwrite destination if it exists

    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Verify source directory exists
        if not os.path.exists(source_dir):
            raise FileNotFoundError(f"Source directory '{source_dir}' does not exist")

        # Manage existing destination directory
        if os.path.exists(destination_dir):
            if overwrite:
                shutil.rmtree(destination_dir)
            else:
                raise FileExistsError(
                    f"Destination directory '{destination_dir}' already exists"
                )

        # Perform the directory copy
        shutil.copytree(source_dir, destination_dir)
        return True

    except PermissionError as e:
        raise PermissionError("Permission denied while copying directory") from e
    except Exception as e:
        raise RuntimeError("An error occurred while copying the directory") from e


def main(
    bsa_dir: Path = typer.Argument(..., help="BSA dir"),
    annotation_file: Path = typer.Argument(..., help="Output file"),
    backup_dir: bool = True,
):
    """ "add annotation to bsa result files"""
    if backup_dir:
        # backup original bsa directory
        back_up_dir = bsa_dir.parent / f"{bsa_dir.name}_backup"
        copy_directory(bsa_dir, back_up_dir)

    anno_df = pd.read_table(annotation_file)
    for data_i in tqdm(bsa_dir.glob(f"*/data/*.gz")):
        df = pd.read_csv(data_i)
        add_anno_df = df.merge(anno_df, how="left", on=["Gene"])
        add_anno_df.to_csv(data_i, index=False)


if __name__ == "__main__":
    typer.run(main)
