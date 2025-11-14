from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional

import typer

TRANSCRIPT_KEYWORDS = {
    "mrna",
    "transcript",
    "lincrna",
    "lnc_rna",
    "ncrna",
    "rrna",
    "trna",
    "snrna",
    "snorna",
    "mirna",
    "pirna",
    "scrna",
    "sirna",
    "primary_transcript",
    "antisense_rna",
    "guide_rna",
}


def parse_attributes(field: str) -> Dict[str, str]:
    """Parse the 9th GFF column into a dict of attributes."""
    attributes: Dict[str, str] = {}
    raw_items = field.strip().strip(";").split(";")
    for item in raw_items:
        item = item.strip()
        if not item:
            continue
        key: Optional[str] = None
        value: Optional[str] = None
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
        if key is None or value is None:
            continue
        attributes[key.strip()] = value.strip().strip('"')
    return attributes


def is_transcript_feature(feature_type: str) -> bool:
    feature = feature_type.lower()
    return (
        feature in TRANSCRIPT_KEYWORDS
        or feature.endswith("rna")
        or "transcript" in feature
    )


def main(gff: Path, out: Optional[Path] = typer.Option(None)) -> None:
    target_json = out or gff.with_suffix(".gene-transcripts.json")
    gene_to_transcripts: Dict[str, set[str]] = defaultdict(set)
    known_gene_ids: set[str] = set()

    with gff.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            attributes = parse_attributes(parts[8])

            if feature_type.lower() == "gene":
                gene_id = (
                    attributes.get("gene_id")
                    or attributes.get("ID")
                    or attributes.get("gene")
                )
                if gene_id:
                    known_gene_ids.add(gene_id)
                    gene_to_transcripts.setdefault(gene_id, set())
                continue

            if not is_transcript_feature(feature_type):
                continue

            transcript_id = (
                attributes.get("transcript_id")
                or attributes.get("transcript")
                or attributes.get("ID")
            )
            if not transcript_id:
                continue

            candidate_genes = []
            if "gene_id" in attributes:
                candidate_genes.append(attributes["gene_id"])
            if "gene" in attributes:
                candidate_genes.append(attributes["gene"])
            if "Parent" in attributes:
                candidate_genes.extend(attributes["Parent"].split(","))

            target_gene = None
            for candidate in candidate_genes:
                candidate = candidate.strip()
                if candidate and (candidate in known_gene_ids or target_gene is None):
                    target_gene = candidate
                    if candidate in known_gene_ids:
                        break

            if target_gene:
                gene_to_transcripts[target_gene].add(transcript_id)

    serialized = {
        gene: sorted(transcripts)
        for gene, transcripts in gene_to_transcripts.items()
        if transcripts
    }

    with target_json.open("w") as out_handle:
        json.dump(serialized, out_handle, indent=2, sort_keys=True)

    typer.echo(f"Wrote {len(serialized)} genes to {target_json}")


if __name__ == "__main__":
    typer.run(main)
