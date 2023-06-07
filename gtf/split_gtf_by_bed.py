from HTSeq import GFF_Reader
import sys

if not len(sys.argv) == 3:
    print "python " + sys.argv[0] + "split_bed gtf_file\n"
    sys.exit(0)

bed = sys.argv[1]
gtf = sys.argv[2]

chr_dict = dict()

with open(bed) as bed_inf:
    for eachline in bed_inf:
        eachline_inf = eachline.strip().split()
        chrom = eachline_inf[0]
        start = int(eachline_inf[1]) + 1
        end = int(eachline_inf[2])
        split_chr = eachline_inf[3]
        chr_dict.setdefault(chrom, {})[(start, end)] = split_chr


for eachline in GFF_Reader(gtf):
    chrom = eachline.iv.chrom
    start = eachline.iv.start + 1
    end = eachline.iv.end
    if chrom in chr_dict:
        for each_inter in chr_dict[chrom]:
            if start >= each_inter[0] and end <= each_inter[1]:
                new_chr = chr_dict[chrom][each_inter]
                new_start = start - each_inter[0] + 1
                new_end = end - each_inter[0] + 1
    else:
        new_chr = chrom
        new_start = start
        new_end = end
    output_line = eachline.get_gff_line().strip().split('\t')
    output_line[0] = new_chr
    output_line[3] = str(new_start)
    output_line[4] = str(new_end)
    print '{l}'.format(l='\t'.join(output_line))
