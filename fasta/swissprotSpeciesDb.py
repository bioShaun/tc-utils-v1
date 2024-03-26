import typer
from Bio import SeqIO
from pathlib import Path
import pandas as pd
from typing import Dict, Optional
import re
from tqdm import tqdm
import gzip
from loguru import logger


def parseSwissprotDescription(description: str) -> Optional[Dict[str, str]]:
    pattern = re.compile("(.*) OS=(.*) OX=(\d+) GN=(.*) PE=(\d+) SV=(\d+)")
    match_res = pattern.match(description)
    if match_res:
        anno, os, ox, gn, pe, sv = match_res.groups()
        swissprot_id = anno.split()[0].split("|")[1]
        swissprot_description = " ".join(anno.split()[1:])
        return {
            "swissprot_id": swissprot_id,
            "swissprot_description": f"{swissprot_description} ({os})",
            "swissprot_taxid": ox,
            "swissprot_name": gn,
        }


def main(uniref_fasta: Path, tax_table: Path, out_prefix: Path) -> None:
    tax_df = pd.read_table(tax_table, index_col=0)
    species_fa_list = []
    fa_info_list = []

    with gzip.open(uniref_fasta, "rt") as handle:
        for record in tqdm(SeqIO.parse(handle, "fasta")):
            record_info = parseSwissprotDescription(record.description)
            if record_info:
                tax_id = int(record_info["ox"])
                if tax_id in tax_df.index:
                    species_fa_list.append(record)
                    fa_info_list.append(record_info)

    SeqIO.write(species_fa_list, out_prefix.with_suffix(".fa"), "fasta")
    if fa_info_list:
        fa_info_df = pd.DataFrame(fa_info_list)
        fa_info_df.to_csv(
            out_prefix.with_suffix(".tsv"),
            index=False,
            sep="\t",
        )
    else:
        logger.warning("No uniref species found")


if __name__ == "__main__":
    typer.run(main)
