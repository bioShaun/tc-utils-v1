#!/usr/bin/env python3
import typer
from cyvcf2 import VCF
from pyfaidx import Fasta

app = typer.Typer(help="Generate probe FASTA from VCF and genome.")


@app.command()
def generate_probes(
    vcf_file: str = typer.Argument(..., help="Input VCF file"),
    genome_fasta: str = typer.Argument(..., help="Reference genome FASTA"),
    probe_size: int = typer.Argument(..., help="Probe size (bp)"),
    output_fasta: str = typer.Argument(..., help="Output probe FASTA file"),
):
    half_size = probe_size // 2
    genome = Fasta(genome_fasta)

    with open(output_fasta, "w") as out_fasta:
        for variant in VCF(vcf_file):
            chrom = variant.CHROM
            pos = variant.POS  # 1-based
            start = max(pos - half_size, 1)
            end = pos + half_size
            try:
                seq = genome[chrom][
                    start - 1 : end
                ].seq  # pyfaidx uses 0-based indexing
                header = f">{chrom}_{pos}"
                out_fasta.write(f"{header}\n{seq}\n")
            except KeyError:
                typer.echo(
                    f"Warning: chromosome {chrom} not found in genome.", err=True
                )
            except Exception as e:
                typer.echo(f"Error processing {chrom}:{pos} -> {e}", err=True)


if __name__ == "__main__":
    app()
