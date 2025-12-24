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
    
    @given(
        ref=st.text(min_size=1, max_size=1, alphabet='ATCG'),
        alt=st.text(min_size=1, max_size=1, alphabet='ATCG')
    )
    @settings(max_examples=50)
    def test_variant_type_classification_accuracy_snp(self, ref, alt):
        """Property 18: Variant Type Classification Accuracy - SNP cases.
        
        Feature: vcf-processor-optimization, Property 18: Variant Type Classification Accuracy
        Validates: Requirements 8.2, 8.3, 8.4, 8.5, 8.6
        """
        assume(ref != alt)  # Ensure it's not a REF type
        
        # Single nucleotide substitution should be classified as SNP
        variant_type = VariantTransformer.classify_variant_type(ref, alt)
        assert variant_type == "SNP"
    
    @given(
        ref=st.text(min_size=1, max_size=10, alphabet='ATCG'),
        alt=st.text(min_size=1, max_size=10, alphabet='ATCG')
    )
    @settings(max_examples=50)
    def test_variant_type_classification_accuracy_indel(self, ref, alt):
        """Property 18: Variant Type Classification Accuracy - INDEL cases.
        
        Feature: vcf-processor-optimization, Property 18: Variant Type Classification Accuracy
        Validates: Requirements 8.2, 8.3, 8.4, 8.5, 8.6
        """
        assume(len(ref) != len(alt))  # Different lengths = INDEL
        assume(ref != alt)  # Ensure it's not a REF type
        
        # Different length variants should be classified as INDEL
        variant_type = VariantTransformer.classify_variant_type(ref, alt)
        assert variant_type == "INDEL"
    
    @given(
        length=st.integers(min_value=2, max_value=5),
        ref=st.text(min_size=2, max_size=5, alphabet='ATCG'),
        alt=st.text(min_size=2, max_size=5, alphabet='ATCG')
    )
    @settings(max_examples=50)
    def test_variant_type_classification_accuracy_mnp(self, length, ref, alt):
        """Property 18: Variant Type Classification Accuracy - MNP cases.
        
        Feature: vcf-processor-optimization, Property 18: Variant Type Classification Accuracy
        Validates: Requirements 8.2, 8.3, 8.4, 8.5, 8.6
        """
        # Ensure same length and different sequences
        assume(len(ref) == len(alt) and len(ref) > 1)
        assume(ref != alt)  # Ensure it's not a REF type
        
        # Same length multi-nucleotide variants should be classified as MNP
        variant_type = VariantTransformer.classify_variant_type(ref, alt)
        assert variant_type == "MNP"
    
    @given(
        ref=st.text(min_size=1, max_size=10, alphabet='ATCG')
    )
    @settings(max_examples=30)
    def test_variant_type_classification_accuracy_ref(self, ref):
        """Property 18: Variant Type Classification Accuracy - REF cases.
        
        Feature: vcf-processor-optimization, Property 18: Variant Type Classification Accuracy
        Validates: Requirements 8.2, 8.3, 8.4, 8.5, 8.6
        """
        # Identical ref and alt should be classified as REF
        variant_type = VariantTransformer.classify_variant_type(ref, ref)
        assert variant_type == "REF"
    
    def test_variant_type_classification_accuracy_symbolic(self):
        """Property 18: Variant Type Classification Accuracy - Symbolic alleles.
        
        Feature: vcf-processor-optimization, Property 18: Variant Type Classification Accuracy
        Validates: Requirements 8.2, 8.3, 8.4, 8.5, 8.6
        """
        # Symbolic deletion should be classified as INDEL
        variant_type = VariantTransformer.classify_variant_type("A", "*")
        assert variant_type == "INDEL"
    
    @given(
        ref=st.text(min_size=1, max_size=5, alphabet='ATCG'),
        alts=st.lists(
            st.text(min_size=1, max_size=5, alphabet='ATCG'),
            min_size=2, max_size=4
        )
    )
    @settings(max_examples=50)
    def test_multi_allelic_type_combination_sorting(self, ref, alts):
        """Property 19: Multi-Allelic Type Combination - Sorting order.
        
        Feature: vcf-processor-optimization, Property 19: Multi-Allelic Type Combination
        Validates: Requirements 8.7
        """
        # Get combined type
        combined_type = VariantTransformer.classify_multi_allelic_types(ref, alts)
        
        # Should be a string
        assert isinstance(combined_type, str)
        
        # If multiple types, should be pipe-separated
        if "|" in combined_type:
            types = combined_type.split("|")
            
            # Should be in the correct order: SNP, INDEL, MNP, REF
            type_order = ["SNP", "INDEL", "MNP", "REF"]
            
            # Check that types appear in the correct order
            last_index = -1
            for variant_type in types:
                current_index = type_order.index(variant_type)
                assert current_index > last_index, f"Types not in correct order: {combined_type}"
                last_index = current_index
        
        # Verify each component type is valid
        valid_types = {"SNP", "INDEL", "MNP", "REF", "UNKNOWN"}
        for variant_type in combined_type.split("|"):
            assert variant_type in valid_types
    
    def test_multi_allelic_type_combination_specific_cases(self):
        """Property 19: Multi-Allelic Type Combination - Specific test cases.
        
        Feature: vcf-processor-optimization, Property 19: Multi-Allelic Type Combination
        Validates: Requirements 8.7
        """
        # Test case with SNP and INDEL
        combined = VariantTransformer.classify_multi_allelic_types("A", ["T", "AT"])
        assert combined == "SNP|INDEL"
        
        # Test case with all types - need to include a real MNP (same length >1)
        combined = VariantTransformer.classify_multi_allelic_types("AT", ["GC", "A", "ATG", "AT"])
        assert combined == "INDEL|MNP|REF"
        
        # Test case with only one type
        combined = VariantTransformer.classify_multi_allelic_types("A", ["T", "G"])
        assert combined == "SNP"
        
        # Test case with REF only
        combined = VariantTransformer.classify_multi_allelic_types("A", ["A"])
        assert combined == "REF"
        
        # Test case with empty alts
        combined = VariantTransformer.classify_multi_allelic_types("A", [])
        assert combined == "REF"
    
    @given(
        ref=st.text(min_size=1, max_size=3, alphabet='ATCG'),
        alts=st.lists(
            st.one_of(
                st.text(min_size=1, max_size=1, alphabet='ATCG'),  # SNPs
                st.text(min_size=2, max_size=5, alphabet='ATCG'),  # Potential INDELs/MNPs
                st.just("*")  # Symbolic deletion
            ),
            min_size=1, max_size=3
        )
    )
    @settings(max_examples=30)
    def test_multi_allelic_type_combination_consistency(self, ref, alts):
        """Property 19: Multi-Allelic Type Combination - Consistency check.
        
        Feature: vcf-processor-optimization, Property 19: Multi-Allelic Type Combination
        Validates: Requirements 8.7
        """
        # Get individual classifications
        individual_types = set()
        for alt in alts:
            variant_type = VariantTransformer.classify_variant_type(ref, alt)
            if variant_type != "UNKNOWN":
                individual_types.add(variant_type)
        
        # Get combined classification
        combined_type = VariantTransformer.classify_multi_allelic_types(ref, alts)
        
        if not individual_types:
            # If no valid types found, should return UNKNOWN
            assert combined_type == "UNKNOWN"
        else:
            # Combined type should contain all individual types
            combined_set = set(combined_type.split("|"))
            assert combined_set == individual_types


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