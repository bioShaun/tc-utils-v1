import typer
from pathlib import Path


def generate_fq_file_map(data_dir: Path, ) -> None:
    pass


def main(
    data_path: Path,
    sample_file: Path,
    out_dir: Path,
    generate: bool = typer.Option(False, "--generate"),
    run: bool = typer.Option(False, "--run"),
    mode: str = typer.Option("link", "--mode"),
) -> None:
    pass


if __name__ == "__main__":
    typer.run(main)
