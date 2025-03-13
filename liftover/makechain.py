from pathlib import Path

import delegator
import typer


def make_chain(ref: Path, query: Path, outdir: Path, threads: int = 16) -> None:
    """Make a chain file from a reference and query file."""
    ref_name = ref.stem
    query_name = query.stem
    outdir.mkdir(exist_ok=True)
    out_paf = outdir / f"{ref_name}_to_{query_name}.paf"
    minimap_cmd = f"minimap2 -t {threads} -cx asm5 --cs {query} {ref} > {out_paf}"
    delegator.run(minimap_cmd)
    chain_file = outdir / f"{ref_name}_to_{query_name}.chain"
    paf2chain_cmd = f"transanno minimap2chain {out_paf} --output {chain_file}"
    delegator.run(paf2chain_cmd)


if __name__ == "__main__":
    typer.run(make_chain)
