import gzip
from enum import Enum
from pathlib import Path
from typing import TextIO

import typer


class VaType(str, Enum):
    SNP = "snp"
    INDEL = "indel"


def filter_no_variant_vcf(
    input_vcf: Path, output_vcf: Path, variant_type: VaType = VaType.SNP
) -> None:
    """
    Filter a VCF file and write the filtered results to an output VCF file.
    """

    with open_vcf(input_vcf) as input_vcf_file, output_vcf.open("w") as output_vcf_file:
        for line in input_vcf_file:
            if line.startswith("#"):
                output_vcf_file.write(line)
                continue

            columns = line.split("\t")
            if columns[4] == ".":
                if variant_type == VaType.SNP and len(columns[3]) > 1:
                    continue
                if variant_type == VaType.INDEL and len(columns[3]) == 1:
                    continue

            output_vcf_file.write(line)


def open_vcf(vcf_path: Path) -> TextIO:
    """
    Open a VCF file. If the file is compressed with gzip, use the gzip module to open it.
    """

    if vcf_path.suffix == ".gz":
        return gzip.open(vcf_path, "rt")
    else:
        return vcf_path.open("r")


if __name__ == "__main__":
    typer.run(filter_no_variant_vcf)
