import gzip
from pathlib import Path

import pytest

from panel.filter_vcf_cluster import filter_stream

# VCF 头部常量
HEADER = """
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"""


def write_vcf(lines, path: Path, gzipped=False):
    """
    将变异行写入临时 VCF 文件，同时添加头部信息，可选是否压缩。
    """
    mode, writer = ("wt", gzip.open) if gzipped else ("w", open)
    with writer(path, mode) as f:
        f.write(HEADER)
        f.writelines(lines)


def read_vcf(path: Path):
    """
    读取 VCF 文件内容，自动处理是否压缩，返回行列表（已去除换行符）。
    """
    if path.suffix in [".gz", ".bgz"]:
        with gzip.open(path, "rt") as f:
            return [l.strip() for l in f.readlines()]
    return [l.strip() for l in Path(path).read_text().splitlines()]


def test_no_clusters(tmp_path):
    """
    当所有 SNP 间距大于窗口大小时，所有记录应被保留。
    """
    lines = [
        "1\t10\t.\tA\tG\t.\t.\t.\n",
        "1\t30\t.\tC\tT\t.\t.\t.\n",
        "1\t60\t.\tG\tA\t.\t.\t.\n",
    ]
    infile = tmp_path / "in.vcf"
    outfile = tmp_path / "out.vcf"
    write_vcf(lines, infile, gzipped=False)

    filter_stream(infile, outfile, window_size=10, snp_count=2)
    output = read_vcf(outfile)

    expected = [l for l in HEADER.strip().splitlines()] + [l.strip() for l in lines]
    assert output == expected


def test_simple_cluster(tmp_path):
    """
    三个 SNP 在窗口范围内时应被删除，后续独立 SNP 保留。
    """
    lines = [
        "1\t100\t.\tA\tG\t.\t.\t.\n",
        "1\t105\t.\tC\tT\t.\t.\t.\n",
        "1\t110\t.\tG\tA\t.\t.\t.\n",
        "1\t200\t.\tT\tC\t.\t.\t.\n",
    ]
    infile = tmp_path / "in.vcf.gz"
    outfile = tmp_path / "out.vcf"
    write_vcf(lines, infile, gzipped=True)

    filter_stream(infile, outfile, window_size=10, snp_count=2)
    output = read_vcf(outfile)

    expected = [l for l in HEADER.strip().splitlines()] + ["1\t200\t.\tT\tC\t.\t.\t."]
    assert output == expected


def test_edge_window(tmp_path):
    """
    对于 snp_count=1，缓冲区大小为 2：当两个变异位置差等于窗口大小时，应被删除。
    """
    lines = [
        "1\t100\t.\tA\tG\t.\t.\t.\n",
        "1\t110\t.\tC\tT\t.\t.\t.\n",
        "1\t120\t.\tG\tA\t.\t.\t.\n",
    ]
    infile = tmp_path / "in.vcf"
    outfile = tmp_path / "out.vcf"
    write_vcf(lines, infile, gzipped=False)

    filter_stream(infile, outfile, window_size=10, snp_count=1)
    output = read_vcf(outfile)

    # 所有记录均被删除，头部保留
    expected = HEADER.strip().splitlines()
    assert output == expected
