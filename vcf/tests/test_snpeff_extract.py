from pathlib import Path

import csv

from vcf.snpeff_extract import extract_best_annotations


def _write_vcf(tmp_path: Path) -> Path:
    vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=ANN,Number=.,Type=String,Description="ANN">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\t.\tA\tG\t.\tPASS\tANN=G|missense_variant|MODERATE|gene1
1\t200\t.\tC\tT\t.\tPASS\tANN=C|stop_gained|HIGH|gene2,A|synonymous_variant|LOW|gene2
1\t300\t.\tG\tA\t.\tPASS\t.
"""
    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(vcf_content)
    return vcf_path


def test_extracts_highest_priority_annotation(tmp_path: Path) -> None:
    vcf_path = _write_vcf(tmp_path)
    output = tmp_path / "out.tsv"

    extract_best_annotations(vcf_path, output)

    with output.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    assert len(rows) == 2  # record without ANN is skipped
    first, second = rows

    assert first["chrom"] == "1"
    assert first["pos"] == "100"
    assert first["effect"] == "missense_variant"
    assert first["impact"] == "MODERATE"

    assert second["chrom"] == "1"
    assert second["pos"] == "200"
    # HIGH is preferred over LOW even if LOW appears first in list
    assert second["effect"] == "stop_gained"
    assert second["impact"] == "HIGH"
