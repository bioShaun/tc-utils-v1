"""Pytest configuration and fixtures for VCF processor tests."""

import tempfile
from pathlib import Path
from typing import Generator

import pytest


@pytest.fixture
def temp_dir() -> Generator[Path, None, None]:
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield Path(tmp_dir)


@pytest.fixture
def sample_vcf_content() -> str:
    """Sample VCF content for testing."""
    return """##fileformat=VCFv4.2
##contig=<ID=chr1,length=249250621>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3
chr1	100	.	A	T	60	PASS	.	GT	0/0	0/1	1/1
chr1	200	.	G	C	60	PASS	.	GT	0/1	1/1	0/0
chr1	300	.	AT	A	60	PASS	.	GT	0/0	0/1	./.
chr1	400	.	C	CG	60	PASS	.	GT	1/1	0/0	0/1
"""


@pytest.fixture
def sample_target_ids() -> str:
    """Sample target IDs for testing."""
    return """chr1_100
chr1_200
chr1_300
"""


@pytest.fixture
def sample_vcf_file(temp_dir: Path, sample_vcf_content: str) -> Path:
    """Create a sample VCF file for testing."""
    vcf_file = temp_dir / "sample.vcf"
    vcf_file.write_text(sample_vcf_content)
    return vcf_file


@pytest.fixture
def sample_target_file(temp_dir: Path, sample_target_ids: str) -> Path:
    """Create a sample target IDs file for testing."""
    target_file = temp_dir / "targets.txt"
    target_file.write_text(sample_target_ids)
    return target_file