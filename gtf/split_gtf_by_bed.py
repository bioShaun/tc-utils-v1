#!/usr/bin/env python3
import sys
from HTSeq import GFF_Reader


def load_intervals(bed_path):
    chr_dict = {}
    with open(bed_path) as bed_inf:
        for eachline in bed_inf:
            chrom, start, end, split_chr = eachline.strip().split()[:4]
            chr_dict.setdefault(chrom, {})[(int(start) + 1, int(end))] = split_chr
    return chr_dict


def main(bed, gtf):
    chr_dict = load_intervals(bed)

    for eachline in GFF_Reader(gtf):
        chrom = eachline.iv.chrom
        start = eachline.iv.start + 1
        end = eachline.iv.end

        new_chr, new_start, new_end = chrom, start, end
        if chrom in chr_dict:
            for each_inter in chr_dict[chrom]:
                if start >= each_inter[0] and end <= each_inter[1]:
                    new_chr = chr_dict[chrom][each_inter]
                    new_start = start - each_inter[0] + 1
                    new_end = end - each_inter[0] + 1
                    break

        output_line = eachline.get_gff_line().strip().split("\t")
        output_line[0] = new_chr
        output_line[3] = str(new_start)
        output_line[4] = str(new_end)
        print("\t".join(output_line))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"python {sys.argv[0]} split_bed gtf_file\n")
        sys.exit(0)

    main(sys.argv[1], sys.argv[2])
