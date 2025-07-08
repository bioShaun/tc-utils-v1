from math import e
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
            if not eachline:
                continue
            if eachline.startswith("#"):
                print(eachline)
                continue
            output_line = eachline.split("\t")
            try:
                chrom, *_, start, end = output_line[:5]
            except ValueError:
                raise ValueError(eachline)
            output_line_cp = output_line[:]
            if chrom in chr_dict:
                for each_inter in chr_dict[chrom]:
                    if int(start) >= each_inter[0]:
                        output_line_cp[0] = chr_dict[chrom][each_inter]
                        output_line_cp[3] = str(int(start) - each_inter[0])
                        if int(end) <= each_inter[1]:
                            output_line_cp[4] = str(int(end) - each_inter[0])
                        else:
                            if int(each_inter[0]) != 0:
                                raise ValueError(f"超过染色体长度, {eachline}")
                            output_line_cp[4] = str(each_inter[1])
                        output_str = "\t".join(output_line_cp)
                        print(output_str)
                    else:
                        if int(end) > each_inter[0]:
                            output_line_cp2 = output_line[:]
                            output_line_cp2[0] = chr_dict[chrom][each_inter]
                            output_line_cp2[3] = str(0)
                            output_line_cp2[4] = str(int(end) - each_inter[0])
                            output_str = "\t".join(output_line_cp2)
                            print(output_str)


if __name__ == "__main__":
    typer.run(main)
