from pathlib import Path

import typer


def main(bed: Path, gtf: Path):
    chr_dict = dict()

    with open(bed) as bed_inf:
        for eachline in bed_inf:
            eachline_inf = eachline.strip().split()
            chrom = eachline_inf[0]
            start = int(eachline_inf[1])
            end = int(eachline_inf[2])
            split_chr = eachline_inf[3]
            chr_dict.setdefault(chrom, {})[(start, end)] = split_chr

    with gtf.open() as gtf_info:
        for eachline in gtf_info:
            eachline = eachline.strip()
            if eachline.startswith("#"):
                print(eachline)
                continue
            output_line = eachline.split("\t")
            chrom, *_, start, end = output_line[:5]
            if chrom in chr_dict:
                for each_inter in chr_dict[chrom]:
                    if start >= each_inter[0] and end <= each_inter[1]:
                        output_line[0] = chr_dict[chrom][each_inter]
                        output_line[3] = int(start) - each_inter[0]
                        output_line[4] = int(end) - each_inter[0]
            output_str = "\t".join(output_line)
            print(output_str)


if __name__ == "__main__":
    typer.run(main)
