import pysam

samfile = pysam.AlignmentFile("aln.sam", "r")

for read in samfile.fetch(until_eof=True):
    # 统计CIGAR中M的长度
    cigar_len = sum(length for (op, length) in read.cigartuples if op == 0)
    # 获取XM（mismatch 数）
    xm = read.get_tag("XM") if read.has_tag("XM") else 0
    effective_len = cigar_len - xm
    print(read.query_name, effective_len)
