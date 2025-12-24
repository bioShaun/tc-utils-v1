"""Property-based tests for variant filtering functionality."""

import tempfile
from pathlib import Path
from typing import Set

import pytest
from hypothesis import given, strategies as st, assume, settings, HealthCheck
from hypothesis.stateful import RuleBasedStateMachine, rule, initialize, invariant

from ..variant_filter import VariantFilter, load_target_ids, validate_target_file
from ..config import VariantInfo


# Helper strategies
def valid_variant_id_strategy():
    """Generate valid variant IDs in CHROM_POS format."""
    return st.builds(
        lambda chrom, pos: f"{chrom}_{pos}",
        chrom=st.text(min_size=1, max_size=10, alphabet='chr0123456789XY'),
        pos=st.integers(min_value=1, max_value=1000000)
    )


class TestVariantFilterProperties:
    """Property-based tests for VariantFilter class."""
    
    @given(
        valid_ids=st.sets(valid_variant_id_strategy(), min_size=1, max_size=20),
        invalid_ids=st.sets(
            st.one_of(
                st.just("no_underscore"),
                st.just("_"),
                st.just("chr_"),
                st.just("_123"),
                st.just("chr_0"),
                st.just("chr_-1"),
                st.just("chr_abc"),
            ),
            max_size=5
        )
    )
    @settings(max_examples=20, suppress_health_check=[HealthCheck.filter_too_much])
    def test_invalid_target_id_handling(self, valid_ids: Set[str], invalid_ids: Set[str]):
        """Property 11: Invalid Target ID Handling.
        
        Validates: Requirements 6.2
        The system should properly identify and handle invalid target IDs.
        """
        assume(not valid_ids.intersection(invalid_ids))  # Ensure no overlap
        
        all_ids = valid_ids.union(invalid_ids)
        
        # Create filter with mixed valid/invalid IDs
        variant_filter = VariantFilter(all_ids)
        
        # Validate target IDs
        detected_invalid = set(variant_filter.validate_target_ids())
        
        # All invalid IDs should be detected
        assert invalid_ids.issubset(detected_invalid)
        
        # No valid IDs should be marked as invalid
        assert not valid_ids.intersection(detected_invalid)
    
    @given(target_ids=st.sets(valid_variant_id_strategy(), min_size=1, max_size=10))
    @settings(max_examples=10, suppress_health_check=[HealthCheck.filter_too_much])
    def test_target_id_file_loading_consistency(self, target_ids: Set[str]):
        """Test that target ID file loading is consistent and handles edge cases."""
        # Create target ID file with comments and empty lines
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("# This is a comment\n")
            f.write("\n")  # Empty line
            for target_id in target_ids:
                f.write(f"{target_id}\n")
            f.write("\n")  # Trailing empty line
            f.write("# Another comment\n")
            target_file = Path(f.name)
        
        try:
            # Load IDs using function
            loaded_ids = load_target_ids(target_file)
            
            # Should match exactly (comments and empty lines ignored)
            assert loaded_ids == target_ids
            
            # Create filter from file
            variant_filter = VariantFilter.from_file(target_file)
            assert variant_filter.target_ids == target_ids
            
            # Validate file
            assert validate_target_file(target_file) == True
            
        finally:
            target_file.unlink()
    
    @given(
        variants=st.lists(
            st.builds(
                VariantInfo,
                chrom=st.sampled_from(['chr1', 'chr2', 'chr3']),
                pos=st.integers(min_value=100, max_value=1000),
                ref=st.sampled_from(['A', 'T', 'G', 'C']),
                alt=st.lists(st.sampled_from(['A', 'T', 'G', 'C']), min_size=1, max_size=2),
                variant_id=valid_variant_id_strategy(),
                genotypes=st.lists(st.sampled_from(['0/0', '0/1', '1/1', './.', '1/0']), min_size=1, max_size=3)
            ),
            min_size=1,
            max_size=10
        )
    )
    @settings(max_examples=10, suppress_health_check=[HealthCheck.filter_too_much])
    def test_filtering_consistency(self, variants):
        """Test that filtering behavior is consistent and predictable."""
        # Extract all variant IDs
        all_variant_ids = {v.variant_id for v in variants}
        
        # Test with no filter (should pass all)
        no_filter = VariantFilter()
        filtered_all = no_filter.filter_variants(variants)
        assert len(filtered_all) == len(variants)
        assert no_filter.get_stats()['passed_filter'] == len(variants)
        
        # Test with subset filter
        if len(all_variant_ids) > 1:
            subset_ids = set(list(all_variant_ids)[:len(all_variant_ids)//2])
            subset_filter = VariantFilter(subset_ids)
            filtered_subset = subset_filter.filter_variants(variants)
            
            # Should only contain variants with IDs in subset
            filtered_ids = {v.variant_id for v in filtered_subset}
            assert filtered_ids.issubset(subset_ids)
            
            # Count should match
            expected_count = sum(1 for v in variants if v.variant_id in subset_ids)
            assert len(filtered_subset) == expected_count
        
        # Test with non-existent IDs (should pass none)
        nonexistent_filter = VariantFilter({'nonexistent_123', 'fake_456'})
        filtered_none = nonexistent_filter.filter_variants(variants)
        assert len(filtered_none) == 0


class VariantFilterStateMachine(RuleBasedStateMachine):
    """Stateful testing for VariantFilter operations."""
    
    def __init__(self):
        super().__init__()
        self.filter = VariantFilter()
        self.added_ids = set()
        self.test_variants = []
    
    @initialize()
    def setup_test_data(self):
        """Initialize with test variants."""
        self.test_variants = [
            VariantInfo(
                chrom='chr1', pos=100, ref='A', alt=['T'], 
                variant_id='chr1_100', genotypes=['0/1', '1/1']
            ),
            VariantInfo(
                chrom='chr1', pos=200, ref='G', alt=['C'], 
                variant_id='chr1_200', genotypes=['0/0', '0/1']
            ),
            VariantInfo(
                chrom='chr2', pos=150, ref='T', alt=['A'], 
                variant_id='chr2_150', genotypes=['1/1', '0/0']
            ),
        ]
    
    @rule(variant_id=st.sampled_from(['chr1_100', 'chr1_200', 'chr2_150', 'chr3_300']))
    def add_target_id(self, variant_id: str):
        """Add a target ID to the filter."""
        self.filter.add_target_id(variant_id)
        self.added_ids.add(variant_id)
    
    @rule(variant_id=st.sampled_from(['chr1_100', 'chr1_200', 'chr2_150', 'chr3_300']))
    def remove_target_id(self, variant_id: str):
        """Remove a target ID from the filter."""
        self.filter.remove_target_id(variant_id)
        self.added_ids.discard(variant_id)
    
    @rule()
    def clear_target_ids(self):
        """Clear all target IDs."""
        self.filter.clear_target_ids()
        self.added_ids.clear()
    
    @rule()
    def filter_test_variants(self):
        """Filter the test variants and verify results."""
        filtered = self.filter.filter_variants(self.test_variants)
        
        # If no target IDs, should pass all variants
        if not self.added_ids:
            assert len(filtered) == len(self.test_variants)
        else:
            # Should only contain variants with IDs in target set
            filtered_ids = {v.variant_id for v in filtered}
            expected_ids = {v.variant_id for v in self.test_variants if v.variant_id in self.added_ids}
            assert filtered_ids == expected_ids
    
    @invariant()
    def target_count_consistent(self):
        """Target count should match added IDs."""
        assert self.filter.get_target_count() == len(self.added_ids)
    
    @invariant()
    def target_ids_consistent(self):
        """Target IDs should match what we've added."""
        assert self.filter.target_ids == self.added_ids


# Run the stateful test
TestVariantFilterState = VariantFilterStateMachine.TestCase