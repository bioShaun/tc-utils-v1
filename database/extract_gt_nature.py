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
                if gt[0] == 1 and gt[1] == 1:
                    continue
                if gt[0] == -1:
                    ref_allele = ref_allele = "."
                else:
                    ref_allele, alt_allele = gt[:2]
                sample_id = sample_names[n]
                gt_inf.write(
                    f"{variant.CHROM}\t{variant.POS}\t{variant.REF}\t{variant.ALT[0]}\t{sample_id}\t{ref_allele}/{alt_allele}\n"
                )


if __name__ == "__main__":
    typer.run(main)
