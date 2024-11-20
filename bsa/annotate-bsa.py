import pandas as pd
import typer


def main(
    bsa_dir: str = typer.Argument(..., help="BSA dir"),
    annotation_file: str = typer.Argument(..., help="Output file"),
):
    """ "add annotation to bsa result files"""
    pass


if __name__ == "__main__":
    typer.run(main)
