import tempfile

import pysam
import pytest

from draft.eff_len_from_sam import (
    compute_effective_lengths,
    parse_md_mismatches,
    process_sam,
)

# 测试 SAM 内容
SAM_CONTENT = """\
@SQ SN:chr1 LN:1000
read1\t0\tchr1\t1\t255\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tXM:i:2\tMD:Z:8A1
read2\t0\tchr1\t11\t255\t5M5M\t*\t0\t0\tCCCCCCCCCC\tIIIIIIIIII\tXM:i:1\tMD:Z:4G5
read3\t0\tchr1\t21\t255\t6M4M\t*\t0\t0\tTTTTTTTTTT\tIIIIIIIIII\tMD:Z:9C
read4\t0\tchr1\t31\t255\t10M\t*\t0\t0\tGGGGGGGGGG\tIIIIIIIIII
read5\t4\t*\t0\t0\t*\t*\t0\t0\tNNNNNNNNNN\tIIIIIIIIII
"""

@pytest.fixture
def tmp_sam():
    with tempfile.NamedTemporaryFile(mode="w+", suffix=".sam") as f:
        f.write(SAM_CONTENT)
        f.flush()
        yield f.name

def test_parse_md_mismatches():
    assert parse_md_mismatches("43G7C58") == 2
    assert parse_md_mismatches("10") == 0
    assert parse_md_mismatches("5A3T2G") == 3

def test_compute_effective_lengths(tmp_sam):
    samfile = pysam.AlignmentFile(tmp_sam, "r")
    reads = list(samfile.fetch(until_eof=True))

    # read1: 10M, XM=2 -> eff_xm=8, MD mismatch=2 -> eff_md=8
    assert compute_effective_lengths(reads[0]) == (8, 8)
    # read2: 10M, XM=1 -> eff_xm=9, MD mismatch=2 -> eff_md=8
    assert compute_effective_lengths(reads[1]) == (9, 8)
    # read3: no XM, MD mismatch=1 -> eff_md=9
    assert compute_effective_lengths(reads[2]) == (None, 9)
    # read4: no XM no MD -> (None, None)
    assert compute_effective_lengths(reads[3]) == (None, None)
    # read5: unmapped -> (None, None)
    assert compute_effective_lengths(reads[4]) == (None, None)

def test_process_sam(tmp_sam):
    with tempfile.NamedTemporaryFile(mode="w+", suffix=".tsv") as f:
        process_sam(tmp_sam, f.name)
        f.seek(0)
        lines = f.read().splitlines()

    assert lines[0] == "readID\teff_len_XM\teff_len_MD"
    assert lines[1].startswith("read1")
    assert "NA" in lines[-1]  # unmapped read 输出 NA
        f.seek(0)
        lines = f.read().splitlines()

    assert lines[0] == "readID\teff_len_XM\teff_len_MD"
    assert lines[1].startswith("read1")
    assert "NA" in lines[-1]  # unmapped read 输出 NA
