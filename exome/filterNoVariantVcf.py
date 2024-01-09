import gzip
import typer
from pathlib import Path


def main(vcf: Path, out_vcf: Path, va_type: str = "snp") -> None:
    out_vcf_inf = out_vcf.open("w")
    with gzip.open(vcf, "rt") as vcf_info:
        for eachline in vcf_info:
            if eachline.startswith("#"):
                out_vcf_inf.write(eachline)
                continue
            eachline_info = eachline.split("\t")
            if eachline_info[4] == ".":
                if va_type == "snp":
                    if len(eachline_info[3]) > 1:
                        continue
                else:
                    if len(eachline_info[3]) == 1:
                        continue
            out_vcf_inf.write(eachline)


if __name__ == "__main__":
    typer.run(main)
