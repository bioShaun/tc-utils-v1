"""Property-based tests for variant transformation functionality."""

import pytest
from hypothesis import given, strategies as st, assume, settings, HealthCheck
from hypothesis.stateful import RuleBasedStateMachine, rule, initialize, invariant

from ..variant_transformer import VariantTransformer, validate_variant_alleles
from ..config import VariantInfo


# Helper strategies
def valid_nucleotide_sequence():
    """Generate valid nucleotide sequences."""
    return st.text(min_size=1, max_size=20, alphabet='ATCG')


def valid_variant_pair():
    """Generate valid ref/alt pairs for testing."""
    return st.one_of(
        # SNPs
        st.tuples(
            st.text(min_size=1, max_size=1, alphabet='ATCG'),
            st.text(min_size=1, max_size=1, alphabet='ATCG')
        ).filter(lambda x: x[0] != x[1]),
        
        # Insertions (ref shorter than alt)
        st.tuples(
            st.text(min_size=1, max_size=5, alphabet='ATCG'),
            st.text(min_size=2, max_size=10, alphabet='ATCG')
        ).filter(lambda x: len(x[0]) < len(x[1])),
        
        # Deletions (ref longer than alt)
        st.tuples(
            st.text(min_size=2, max_size=10, alphabet='ATCG'),
            st.text(min_size=1, max_size=5, alphabet='ATCG')
        ).filter(lambda x: len(x[0]) > len(x[1])),
        
        # MNPs (same length, different sequence)
        st.tuples(
            st.text(min_size=2, max_size=5, alphabet='ATCG'),
            st.text(min_size=2, max_size=5, alphabet='ATCG')
        ).filter(lambda x: len(x[0]) == len(x[1]) and x[0] != x[1])
    )


class TestVariantTransformerProperties:
    """Property-based tests for VariantTransformer class."""
    
    @given(ref_alt=valid_variant_pair())
    @settings(max_examples=50, suppress_health_check=[HealthCheck.filter_too_much])
    def test_variant_transformation_preservation(self, ref_alt):
        """Property 17: Variant Transformation Preservation.
        
        Validates: Requirements 7.5
        Transformation should preserve the essential variant information.
        """
        ref, alt = ref_alt
        
        transformer = VariantTransformer(strict_mode=False)
        
        # Transform the variant
        result = transformer.transform_one_alt(ref, alt)
        
        # Result should not be empty
        assert result
        assert isinstance(result, str)
        
        # For simple cases, verify expected patterns
        if len(ref) == 1 and len(alt) == 1:
            # SNP should return the alt base
            assert result == alt
        elif len(ref) > len(alt):
            # Deletion should start with 'del' or be the alt sequence
            assert result.startswith('del') or result == alt
        elif len(ref) < len(alt):
            # Insertion should start with 'ins' or be the alt sequence
            assert result.startswith('ins') or result == alt
    
    @given(
        ref=valid_nucleotide_sequence(),
        alt_list=st.lists(valid_nucleotide_sequence(), min_size=1, max_size=5)
    )
    @settings(max_examples=30)
    def test_multiple_alt_consistency(self, ref, alt_list):
        """Test that multiple alternative alleles are handled consistently."""
        transformer = VariantTransformer(strict_mode=False)
        
        # Transform each alt individually
        individual_results = [transformer.transform_one_alt(ref, alt) for alt in alt_list]
        
        # Transform all alts together
        batch_results = transformer.transform_alt_alleles(ref, alt_list)
        
        # Results should be identical
        assert individual_results == batch_results
        
        # String method should also match
        alt_string = ",".join(alt_list)
        string_result = transformer.transform_alt_string(ref, alt_string)
        expected_string = ",".join(batch_results)
        assert string_result == expected_string
    
    @given(
        invalid_inputs=st.one_of(
            st.tuples(st.just(""), st.text(min_size=1, max_size=5, alphabet='ATCG')),  # Empty ref
            st.tuples(st.text(min_size=1, max_size=5, alphabet='ATCG'), st.just("")),  # Empty alt
            st.tuples(st.just(""), st.just("")),  # Both empty
        )
    )
    @settings(max_examples=20)
    def test_transformation_error_context(self, invalid_inputs):
        """Property 12: Transformation Error Context.
        
        Validates: Requirements 6.3
        Errors should provide meaningful context about what went wrong.
        """
        ref, alt = invalid_inputs
        
        transformer_strict = VariantTransformer(strict_mode=True)
        transformer_lenient = VariantTransformer(strict_mode=False)
        
        # Strict mode should raise ValueError with meaningful message
        with pytest.raises(ValueError) as exc_info:
            transformer_strict.transform_one_alt(ref, alt)
        
        error_msg = str(exc_info.value)
        assert "empty" in error_msg.lower()
        
        # Lenient mode should handle gracefully
        result = transformer_lenient.transform_one_alt(ref, alt)
        assert result == alt  # Should return original alt
        
        # Error should be recorded in stats
        stats = transformer_lenient.get_stats()
        assert stats['errors'] > 0
    
    @given(
        variants=st.lists(
            st.builds(
                VariantInfo,
                chrom=st.sampled_from(['chr1', 'chr2']),
                pos=st.integers(min_value=100, max_value=1000),
                ref=valid_nucleotide_sequence(),
                alt=st.lists(valid_nucleotide_sequence(), min_size=1, max_size=3),
                variant_id=st.builds(lambda c, p: f"{c}_{p}", 
                                   st.sampled_from(['chr1', 'chr2']),
                                   st.integers(min_value=100, max_value=1000)),
                genotypes=st.lists(st.sampled_from(['0/0', '0/1', '1/1']), min_size=1, max_size=3)
            ),
            min_size=1,
            max_size=10
        )
    )
    @settings(max_examples=20)
    def test_variant_object_transformation(self, variants):
        """Test transformation of VariantInfo objects."""
        transformer = VariantTransformer(strict_mode=False)
        
        for variant in variants:
            transformed = transformer.transform_variant(variant)
            
            # Basic properties should be preserved
            assert transformed.chrom == variant.chrom
            assert transformed.pos == variant.pos
            assert transformed.variant_id == variant.variant_id
            assert transformed.genotypes == variant.genotypes
            
            # Reference should be single base (first base of original)
            assert len(transformed.ref) == 1
            assert transformed.ref == variant.ref[0]
            
            # Alt alleles should be transformed
            assert len(transformed.alt) == len(variant.alt)
            for orig_alt, trans_alt in zip(variant.alt, transformed.alt):
                assert trans_alt  # Should not be empty
    
    def test_validation_function(self):
        """Test the validation helper function."""
        # Valid inputs
        assert validate_variant_alleles("A", ["T"]) == []
        assert validate_variant_alleles("ATG", ["A", "T"]) == []
        assert validate_variant_alleles("C", ["*"]) == []  # Symbolic allele
        
        # Invalid inputs
        issues = validate_variant_alleles("", ["T"])
        assert any("empty" in issue.lower() for issue in issues)
        
        issues = validate_variant_alleles("A", [])
        assert any("no alternative" in issue.lower() for issue in issues)
        
        issues = validate_variant_alleles("X", ["T"])  # Invalid nucleotide
        assert any("invalid characters" in issue.lower() for issue in issues)


class VariantTransformerStateMachine(RuleBasedStateMachine):
    """Stateful testing for VariantTransformer operations."""
    
    def __init__(self):
        super().__init__()
        self.transformer = VariantTransformer(strict_mode=False)
        self.processed_variants = []
    
    @rule(
        ref=st.text(min_size=1, max_size=5, alphabet='ATCG'),
        alt=st.text(min_size=1, max_size=5, alphabet='ATCG')
    )
    def transform_variant(self, ref: str, alt: str):
        """Transform a variant and record it."""
        result = self.transformer.transform_one_alt(ref, alt)
        self.processed_variants.append((ref, alt, result))
    
    @rule()
    def check_stats_consistency(self):
        """Verify that statistics are consistent with processed variants."""
        stats = self.transformer.get_stats()
        total_processed = len(self.processed_variants)
        total_in_stats = sum(stats.values())
        
        # Total in stats should match processed variants
        assert total_in_stats == total_processed
    
    @rule()
    def get_summary(self):
        """Test summary generation."""
        summary = self.transformer.get_summary()
        assert isinstance(summary, str)
        
        if self.processed_variants:
            assert "Total:" in summary
        else:
            assert "No variants processed" in summary
    
    @rule()
    def reset_stats(self):
        """Reset statistics and verify clean state."""
        self.transformer.reset_stats()
        stats = self.transformer.get_stats()
        
        # All stats should be zero
        assert all(count == 0 for count in stats.values())
        
        # Clear our tracking too
        self.processed_variants.clear()
    
    @invariant()
    def stats_non_negative(self):
        """Statistics should never be negative."""
        stats = self.transformer.get_stats()
        assert all(count >= 0 for count in stats.values())


# Run the stateful test
TestVariantTransformerState = VariantTransformerStateMachine.TestCase