import gzip
from pathlib import Path

import typer


def split_chr_size_inf(bed: Path) -> str:
    split_chr_size_list = []
    with open(bed) as bed_inf:
        for row in bed_inf:
            _, start, end, chrom = row.strip().split()
            chrom_length = int(end) - int(start)
            eachline = f"##contig=<ID={chrom},length={chrom_length}>"
            split_chr_size_list.append(eachline)
    return "\n".join(split_chr_size_list)


def main(bed: Path, vcf: Path):
    chr_dict = dict()

    split_chr_size_str = split_chr_size_inf(bed)

    with open(bed) as bed_inf:
        for eachline in bed_inf:
            eachline_inf = eachline.strip().split()
            chrom = eachline_inf[0]
            start = int(eachline_inf[1])
            end = int(eachline_inf[2])
            split_chr = eachline_inf[3]
            chr_dict.setdefault(chrom, {})[(start, end)] = split_chr

    chr_header_flag = 1
    with gzip.open(vcf, "rt") as gtf_info:
        for eachline in gtf_info:
            eachline = eachline.strip()
            if eachline.startswith("##contig="):
                if chr_header_flag:
                    print(split_chr_size_str)
                    chr_header_flag = 0
                continue
            if eachline.startswith("#"):
                print(eachline)
                continue
            output_line = eachline.split("\t")
            chrom, pos = output_line[:2]
            if chrom in chr_dict:
                for each_inter in chr_dict[chrom]:
                    if int(pos) > each_inter[0] and int(pos) <= each_inter[1]:
                        output_line[0] = chr_dict[chrom][each_inter]
                        output_line[1] = str(int(pos) - each_inter[0])
            output_str = "\t".join(output_line)
            print(output_str)


if __name__ == "__main__":
    typer.run(main)
