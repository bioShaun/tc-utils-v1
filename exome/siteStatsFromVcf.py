from pathlib import Path
from typing import List
import vcfpy
from collections import Counter
import typer
from tqdm import tqdm

import gzip


def get_allele_stats(record: vcfpy.Record) -> List[float]:
    gt_list = [each.data.get("GT").replace("|", "/") for each in record.calls]
    gt_count = Counter(gt_list)
    missing_count = gt_count.get("./.", 0)
    valid_count = len(gt_list) - missing_count
    alt_count = gt_count.get("1/1", 0)
    het_count = gt_count.get("0/1", 0)
    missing_rate = round(missing_count / len(gt_list), 3)
    het_rate = round(het_count / valid_count, 3)
    af = (alt_count + het_count * 0.5) / valid_count
    return [missing_rate, het_rate, af]


def main(vcf: Path, site_stats: Path) -> None:
    with gzip.open(site_stats, "wt") as site_inf:
        reader = vcfpy.Reader.from_path(vcf)
        for record in tqdm(reader):
            if record:
                missing_rate, het_rate, af = get_allele_stats(record)
                alt = ",".join([each.value for each in record.ALT])
                site_inf.write(
                    f'{record.CHROM},{record.POS},"{record.REF},{alt}",{missing_rate},{het_rate},{af}\n'
                )


if __name__ == "__main__":
    typer.run(main)
