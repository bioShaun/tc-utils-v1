from pathlib import Path

import pandas as pd
import typer


def main(data_directory: Path, family_file: Path, output_directory: Path) -> None:
    """Add sample information to STRUCTURE output files."""
    output_directory.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(
        family_file,
        sep=r"\s",
        header=None,
        usecols=[1],
        names=["sample_id"],
    )
    for population_file in data_directory.glob("*.Q"):
        populations = pd.read_csv(
            population_file,
            sep=r"\s",
            header=None,
        )
        populations.columns = [f"population_{int(i)+1}" for i in populations.columns]
        combined = pd.concat([samples, populations], axis=1)
        output_file = output_directory / f"{population_file.name}.xlsx"
        combined.to_excel(output_file, index=False)


if __name__ == "__main__":
    typer.run(main)
