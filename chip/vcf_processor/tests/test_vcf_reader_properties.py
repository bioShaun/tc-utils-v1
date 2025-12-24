"""Property-based tests for VCF reading functionality."""

import tempfile
from pathlib import Path
from typing import Set

import pytest
from hypothesis import given, strategies as st, assume, settings
from hypothesis.stateful import RuleBasedStateMachine, rule, initialize, invariant

from ..vcf_reader import VCFReader, load_target_ids
from ..config import VariantInfo


class TestVCFReaderProperties:
    """Property-based tests for VCFReader class."""
    
    @given(
        chrom=st.text(min_size=1, max_size=10, alphabet=st.characters(whitelist_categories=('Lu', 'Ll', 'Nd'))),
        pos=st.integers(min_value=1, max_value=1000000)
    )
    @settings(max_examples=50)
    def test_variant_id_generation_consistency(self, chrom: str, pos: int):
        """Property 1: Variant ID Generation Consistency.
        
        Validates: Requirements 1.2
        The variant ID should always be generated consistently as CHROM_POS format.
        """
        # Create a minimal VCF content for testing
        vcf_content = f"""##fileformat=VCFv4.2
##contig=<ID={chrom}>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
{chrom}	{pos}	.	A	T	60	PASS	.	GT	0/1
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(vcf_content)
            vcf_path = Path(f.name)
        
        try:
            with VCFReader(vcf_path) as reader:
                variants = list(reader.iter_variants())
                
            # Should have exactly one variant
            assert len(variants) == 1
            variant = variants[0]
            
            # Variant ID should follow CHROM_POS format
            expected_id = f"{chrom}_{pos}"
            assert variant.variant_id == expected_id
            
            # Verify other fields match
            assert variant.chrom == chrom
            assert variant.pos == pos
            
        finally:
            vcf_path.unlink()
    
    @given(
        sample_names=st.lists(
            st.text(min_size=1, max_size=20, alphabet=st.characters(whitelist_categories=('Lu', 'Ll', 'Nd', 'Pc'))),
            min_size=1,
            max_size=10,
            unique=True
        )
    )
    @settings(max_examples=30)
    def test_sample_name_extraction_consistency(self, sample_names):
        """Test that sample names are extracted consistently from VCF header."""
        # Create VCF with specified sample names
        sample_header = '\t'.join(sample_names)
        vcf_content = f"""##fileformat=VCFv4.2
##contig=<ID=chr1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{sample_header}
chr1	100	.	A	T	60	PASS	.	GT	{'	'.join(['0/1'] * len(sample_names))}
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(vcf_content)
            vcf_path = Path(f.name)
        
        try:
            with VCFReader(vcf_path) as reader:
                extracted_names = reader.sample_names
                
            # Sample names should match exactly
            assert extracted_names == sample_names
            assert len(extracted_names) == len(sample_names)
            
        finally:
            vcf_path.unlink()
    
    @given(
        target_ids=st.sets(
            st.text(min_size=1, max_size=20, alphabet=st.characters(whitelist_categories=('Lu', 'Ll', 'Nd', 'Pc'))),
            min_size=1,
            max_size=100
        )
    )
    @settings(max_examples=20)
    def test_target_id_filtering_consistency(self, target_ids: Set[str]):
        """Test that target ID filtering works consistently."""
        # Create target ID file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            for target_id in target_ids:
                f.write(f"{target_id}\n")
            target_file = Path(f.name)
        
        try:
            loaded_ids = load_target_ids(target_file)
            
            # Loaded IDs should match input
            assert loaded_ids == target_ids
            assert len(loaded_ids) == len(target_ids)
            
        finally:
            target_file.unlink()


class VCFReaderStateMachine(RuleBasedStateMachine):
    """Stateful testing for VCFReader operations."""
    
    def __init__(self):
        super().__init__()
        self.vcf_path = None
        self.reader = None
        self.expected_variants = []
    
    @initialize()
    def setup_vcf(self):
        """Initialize with a test VCF file."""
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr2>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2
chr1	100	.	A	T	60	PASS	.	GT	0/1	1/1
chr1	200	.	G	C	60	PASS	.	GT	0/0	0/1
chr2	150	.	T	A	60	PASS	.	GT	1/1	0/0
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(vcf_content)
            self.vcf_path = Path(f.name)
        
        # Expected variants for validation
        self.expected_variants = [
            ("chr1", 100, "A", ["T"], "chr1_100"),
            ("chr1", 200, "G", ["C"], "chr1_200"),
            ("chr2", 150, "T", ["A"], "chr2_150"),
        ]
    
    @rule()
    def open_reader(self):
        """Open VCF reader."""
        if self.reader is None:
            self.reader = VCFReader(self.vcf_path)
            self.reader.__enter__()
    
    @rule()
    def close_reader(self):
        """Close VCF reader."""
        if self.reader is not None:
            self.reader.__exit__(None, None, None)
            self.reader = None
    
    @rule()
    def read_all_variants(self):
        """Read all variants and verify consistency."""
        if self.reader is not None:
            variants = list(self.reader.iter_variants())
            
            # Should match expected number
            assert len(variants) == len(self.expected_variants)
            
            # Verify each variant
            for i, (expected_chrom, expected_pos, expected_ref, expected_alt, expected_id) in enumerate(self.expected_variants):
                variant = variants[i]
                assert variant.chrom == expected_chrom
                assert variant.pos == expected_pos
                assert variant.ref == expected_ref
                assert variant.alt == expected_alt
                assert variant.variant_id == expected_id
    
    @rule(target_id=st.sampled_from(["chr1_100", "chr2_150", "nonexistent"]))
    def read_filtered_variants(self, target_id: str):
        """Read variants with filtering."""
        if self.reader is not None:
            target_ids = {target_id}
            variants = list(self.reader.iter_variants(target_ids))
            
            if target_id in ["chr1_100", "chr2_150"]:
                # Should find the variant
                assert len(variants) == 1
                assert variants[0].variant_id == target_id
            else:
                # Should find no variants
                assert len(variants) == 0
    
    @invariant()
    def sample_names_consistent(self):
        """Sample names should always be consistent."""
        if self.reader is not None:
            names = self.reader.sample_names
            assert names == ["sample1", "sample2"]
            assert len(names) == 2
    
    def teardown(self):
        """Clean up resources."""
        if self.reader is not None:
            self.reader.__exit__(None, None, None)
        if self.vcf_path and self.vcf_path.exists():
            self.vcf_path.unlink()


# Run the stateful test
TestVCFReaderState = VCFReaderStateMachine.TestCase