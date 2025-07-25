import gzip
import re
from pathlib import Path
from typing import Dict, Optional

import pandas as pd
import typer
from Bio import SeqIO
from loguru import logger
from tqdm import tqdm


def parseUnirefDescription(description: str) -> Optional[Dict[str, str]]:
    pattern = re.compile(r"(.*) n=(\d+) Tax=(.*) TaxID=(\d+) RepID=(.*)")
    match_res = pattern.match(description)
    if match_res:
        anno, n, tax, taxid, repid = match_res.groups()
        uniref_id = anno.split()[0]
        uniref_description = " ".join(anno.split()[1:])
        return {
            "uniref_id": repid,
            "uniref_description": f"{uniref_description} ({tax})",
            "uniref_taxid": taxid,
        }


def main(uniref_fasta: Path, tax_table: Path, out_prefix: Path) -> None:
    tax_df = pd.read_table(tax_table, index_col=0)
    species_fa_list = []
    fa_info_list = []

    with gzip.open(uniref_fasta, "rt") as handle:
        for record in tqdm(SeqIO.parse(handle, "fasta")):
            record_info = parseUnirefDescription(record.description)
            if record_info:
                tax_id = int(record_info["uniref_taxid"])
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
            columns=["uniref_id", "uniref_description"],
        )
    else:
        logger.warning("No uniref species found")


if __name__ == "__main__":
    typer.run(main)
