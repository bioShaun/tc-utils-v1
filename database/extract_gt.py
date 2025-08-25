from pathlib import Path

import cyvcf2
import typer
from tqdm import tqdm


def main(vcf_file: Path, gt_file: Path) -> None:
    vcf_obj = cyvcf2.VCF(vcf_file)
    sample_names = vcf_obj.samples
    with gt_file.open("w") as gt_inf:
        for variant in tqdm(vcf_obj):
            for n, gt in enumerate(variant.genotypes):
                if gt[0] != -1:
                    ref_count = variant.gt_ref_depths[n]
                    alt_count = variant.gt_alt_depths[n]
                    sample_id = sample_names[n]
                    ref_allele, alt_allele = gt[:2]
                    gt_inf.write(
                        f"{variant.CHROM}\t{variant.POS}\t{variant.REF}\t{variant.ALT[0]}\t.\t{sample_id}={ref_allele}/{alt_allele};{ref_count},{alt_count}\n"
                    )


if __name__ == "__main__":
    typer.run(main)
