from typing import List, Tuple

import typer
from cyvcf2 import VCF
from pyfaidx import Fasta
from tqdm import tqdm

app = typer.Typer(help="Generate probe FASTA from VCF and genome.")


def extract_probe(
    chrom: str, pos: int, genome: Fasta, probe_size: int
) -> Tuple[str, str]:
    """
    Generate probe sequence for a single VCF position.
    Returns (header, sequence)
    """
    half_size = probe_size // 2
    start = max(pos - half_size, 1)
    end = pos + half_size
    seq = genome[chrom][start - 1 : end].seq  # pyfaidx uses 0-based
    header = f">{chrom}_{pos}"
    return header, seq


def generate_probes_from_vcf(
    vcf_file: str, genome_fasta: str, probe_size: int
) -> List[Tuple[str, str]]:
    genome = Fasta(genome_fasta)
    probes = []
    for variant in tqdm(VCF(vcf_file), desc="Processing VCF"):
        try:
            header, seq = extract_probe(variant.CHROM, variant.POS, genome, probe_size)
            probes.append((header, seq))
        except KeyError:
            typer.echo(
                f"Warning: chromosome {variant.CHROM} not found in genome.", err=True
            )
        except Exception as e:
            typer.echo(
                f"Error processing {variant.CHROM}:{variant.POS} -> {e}", err=True
            )
    return probes


def write_fasta(probes: List[Tuple[str, str]], output_fasta: str):
    with open(output_fasta, "w") as out_f:
        for header, seq in probes:
            out_f.write(f"{header}\n{seq}\n")


@app.command()
def main(
    vcf_file: str = typer.Argument(..., help="Input VCF file"),
    genome_fasta: str = typer.Argument(..., help="Reference genome FASTA"),
    probe_size: int = typer.Argument(..., help="Probe size (bp)"),
    output_fasta: str = typer.Argument(..., help="Output probe FASTA file"),
):
    probes = generate_probes_from_vcf(vcf_file, genome_fasta, probe_size)
    write_fasta(probes, output_fasta)


if __name__ == "__main__":
    app()
    app()
