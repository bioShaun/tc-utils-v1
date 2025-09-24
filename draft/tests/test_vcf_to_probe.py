import tempfile

import pytest

from draft.vcf_to_probe import generate_probes_from_vcf, write_fasta

# 可控小型基因组
FASTA_CONTENT = """>chr1
ACGTACGTACGTACGTACGT
"""

# 对应 VCF，pos 为 1-based
VCF_CONTENT = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	3	.	G	A	.	.	.
chr1	8	.	T	C	.	.	.
chr1	1	.	A	T	.	.	.
"""


def create_temp_file(content: str, suffix: str):
    f = tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=suffix)
    f.write(content)
    f.flush()
    return f.name


def test_probe_sequences_match_genome():
    vcf_file = create_temp_file(VCF_CONTENT, ".vcf")
    fasta_file = create_temp_file(FASTA_CONTENT, ".fa")
    output_file = tempfile.NamedTemporaryFile(
        mode="w+", delete=False, suffix=".fa"
    ).name

    probe_size = 5
    probes = generate_probes_from_vcf(vcf_file, fasta_file, probe_size)
    write_fasta(probes, output_file)

    # 从参考基因组中手动计算期望序列
    genome_seq = "ACGTACGTACGTACGTACGT"
    expected_sequences = {
        "chr1_3": genome_seq[0:5],  # pos=3, half=2, start=1,end=5
        "chr1_8": genome_seq[5:10],  # pos=8, start=6,end=10
        "chr1_1": genome_seq[0:3],  # pos=1, start=1,end=3
    }

    with open(output_file) as f:
        lines = [line.strip() for line in f.readlines()]

    # lines[0]=header, lines[1]=seq, lines[2]=header, ...
    for i in range(0, len(lines), 2):
        header = lines[i]
        seq = lines[i + 1]
        assert seq == expected_sequences[header], f"{header} sequence mismatch"
        seq = lines[i + 1]
        assert seq == expected_sequences[header], f"{header} sequence mismatch"
