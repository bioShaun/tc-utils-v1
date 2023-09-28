from pathlib import Path
import typer


def main(gtf: Path) -> None:
    out_gtf = gtf.with_suffix("addGene.gtf")
    with out_gtf.open("w") as out:
        with gtf.open("r") as in_gtf:
            for line in in_gtf:
                if "gene_id" in line:
                    out.write(line)
                else:
                    gene_id = line.strip().split()[:-1]
                    out.write(f"{line.strip()} gene_id {gene_id}")


if __name__ == "__main__":
    typer.run(main)
