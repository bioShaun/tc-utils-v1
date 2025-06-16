#!/usr/bin/env python3

import itertools
from pathlib import Path
from typing import List

import typer

app = typer.Typer(help="Generate pairwise combinations from a list of sample names.")


def read_sample_list(filename: Path) -> List[str]:
    """Read sample names from a file."""
    with open(filename, "r") as f:
        # Strip whitespace and filter out empty lines
        samples = [line.strip() for line in f if line.strip()]
    return samples


def generate_pairs(samples: List[str]) -> List[tuple]:
    """Generate all possible pairwise combinations."""
    return list(itertools.combinations(samples, 2))


@app.command()
def main(
    input_file: Path = typer.Argument(
        ...,
        help="Path to the file containing sample names (one per line)",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
    ),
    output_file: Path = typer.Option(
        None,
        "--output",
        "-o",
        help="Path to save the output (if not specified, prints to stdout)",
        file_okay=True,
        dir_okay=False,
        writable=True,
    ),
):
    """
    Generate all possible pairwise combinations from a list of sample names.

    The input file should contain one sample name per line.
    Output format: sampleA<tab>sampleB
    """
    try:
        samples = read_sample_list(input_file)
        if not samples:
            typer.echo("Error: No samples found in the input file.", err=True)
            raise typer.Exit(code=1)

        pairs = generate_pairs(samples)

        output = f"Found {len(samples)} samples, generating {len(pairs)} pairwise combinations:\n"
        output += "\n".join(f"{pair[0]}\t{pair[1]}" for pair in pairs)

        if output_file:
            output_file.write_text(output)
            typer.echo(f"Results saved to {output_file}")
        else:
            typer.echo(output)

    except Exception as e:
        typer.echo(f"Error: {str(e)}", err=True)
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
