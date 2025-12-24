"""Property-based tests for genotype conversion functionality."""

import pandas as pd
from hypothesis import given, strategies as st, assume, settings, HealthCheck
from hypothesis.stateful import RuleBasedStateMachine, rule, initialize, invariant

from ..genotype_converter import GenotypeConverter, BatchProcessor
from ..config import VariantInfo


# Helper strategies
def valid_genotype_strategy():
    """Generate valid VCF genotype strings."""
    return st.one_of(
        st.just("./."),  # Missing
        st.just("0/0"),  # Homozygous ref
        st.just("0/1"),  # Heterozygous
        st.just("1/0"),  # Heterozygous (reverse)
        st.just("1/1"),  # Homozygous alt
        st.just("0/2"),  # Heterozygous with second alt
        st.just("1/2"),  # Heterozygous between alts
        st.just("2/2"),  # Homozygous second alt
        st.just("0|1"),  # Phased heterozygous
        st.just("1|0"),  # Phased heterozygous (reverse)
    )


def valid_allele_strategy():
    """Generate valid DNA sequences."""
    return st.text(min_size=1, max_size=10, alphabet='ATCG')


class TestGenotypeConverterProperties:
    """Property-based tests for GenotypeConverter class."""
    
    @given(
        genotype=valid_genotype_strategy(),
        ref=valid_allele_strategy(),
        alt_alleles=st.lists(valid_allele_strategy(), min_size=1, max_size=3)
    )
    @settings(max_examples=100, suppress_health_check=[HealthCheck.filter_too_much])
    def test_output_format_equivalence(self, genotype, ref, alt_alleles):
        """Property 2: Output Format Equivalence.
        
        Validates: Requirements 1.5, 7.1, 7.2
        The output format should be consistent and follow expected patterns.
        """
        converter = GenotypeConverter(miss_fmt="NN", gt_sep="")
        
        result = converter.convert_genotype(genotype, ref, alt_alleles)
        
        # Result should not be empty
        assert result
        assert isinstance(result, str)
        
        # Check specific patterns based on genotype
        if genotype == "./.":
            assert result == "NN"  # Missing format
        elif genotype in ["0/0", "0|0"]:
            # Homozygous ref - should be ref sequence or ref+ref
            if len(ref) > 1:
                assert result == ref
            else:
                assert result == ref or result == f"{ref}{ref}"
        elif genotype in ["1/1", "1|1"]:
            # Homozygous alt - should be alt sequence or alt+alt
            alt1 = alt_alleles[0]
            if len(alt1) > 1:
                assert result == alt1
            else:
                assert result == alt1 or result == f"{alt1}{alt1}"
        elif genotype in ["0/1", "1/0", "0|1", "1|0"]:
            # Heterozygous - should contain both alleles
            alt1 = alt_alleles[0]
            if len(ref) > 1 or len(alt1) > 1:
                assert "/" in result
                assert ref in result and alt1 in result
            else:
                # Single base alleles
                assert ref in result and alt1 in result
    
    @given(
        variants=st.lists(
            st.builds(
                VariantInfo,
                chrom=st.sampled_from(['chr1', 'chr2']),
                pos=st.integers(min_value=100, max_value=1000),
                ref=valid_allele_strategy(),
                alt=st.lists(valid_allele_strategy(), min_size=1, max_size=2),
                variant_id=st.builds(lambda c, p: f"{c}_{p}", 
                                   st.sampled_from(['chr1', 'chr2']),
                                   st.integers(min_value=100, max_value=1000)),
                genotypes=st.lists(valid_genotype_strategy(), min_size=1, max_size=5)
            ),
            min_size=1,
            max_size=10
        ),
        sample_names=st.lists(
            st.text(min_size=1, max_size=10, alphabet=st.characters(whitelist_categories=('Lu', 'Ll', 'Nd'))),
            min_size=1,
            max_size=5
        )
    )
    @settings(max_examples=20)
    def test_batch_processing_consistency(self, variants, sample_names):
        """Test that batch processing produces consistent results."""
        # Ensure genotype count matches sample count for each variant
        for variant in variants:
            if len(variant.genotypes) != len(sample_names):
                # Adjust genotypes to match sample count
                if len(variant.genotypes) < len(sample_names):
                    variant.genotypes.extend(["0/0"] * (len(sample_names) - len(variant.genotypes)))
                else:
                    variant.genotypes = variant.genotypes[:len(sample_names)]
        
        converter = GenotypeConverter()
        
        # Convert individually
        individual_results = []
        for variant in variants:
            converted = converter.convert_variant_genotypes(variant, sample_names)
            individual_results.append(converted)
        
        # Convert as batch
        batch_df = converter.convert_batch(variants, sample_names)
        
        # Results should be consistent
        assert len(batch_df) == len(variants)
        
        for i, (variant, individual_result) in enumerate(zip(variants, individual_results)):
            batch_row = batch_df.iloc[i]
            
            # Check location columns
            assert batch_row['CHROM'] == variant.chrom
            assert batch_row['POS'] == variant.pos
            assert batch_row['REF'] == variant.ref
            
            # Check genotype conversions
            for sample_name in sample_names:
                assert batch_row[sample_name] == individual_result[sample_name]
    
    @given(
        batch_size=st.integers(min_value=1, max_value=10),
        variants=st.lists(
            st.builds(
                VariantInfo,
                chrom=st.just('chr1'),
                pos=st.integers(min_value=100, max_value=200),
                ref=st.just('A'),
                alt=st.just(['T']),
                variant_id=st.builds(lambda p: f"chr1_{p}", st.integers(min_value=100, max_value=200)),
                genotypes=st.lists(st.sampled_from(['0/0', '0/1', '1/1']), min_size=2, max_size=2)
            ),
            min_size=1,
            max_size=20
        )
    )
    @settings(max_examples=20)
    def test_batch_size_respect(self, batch_size, variants):
        """Property 5: Batch Size Respect.
        
        Validates: Requirements 4.3
        Batch processing should respect the specified batch size.
        """
        sample_names = ['sample1', 'sample2']
        converter = GenotypeConverter()
        processor = BatchProcessor(batch_size=batch_size)
        
        batches = list(processor.process_variants(iter(variants), sample_names, converter))
        
        # All batches except possibly the last should have batch_size variants
        for i, batch_df in enumerate(batches[:-1]):
            assert len(batch_df) == batch_size
        
        # Last batch should have remaining variants
        if batches:
            last_batch = batches[-1]
            expected_last_size = len(variants) % batch_size
            if expected_last_size == 0:
                expected_last_size = batch_size
            assert len(last_batch) == expected_last_size
        
        # Total variants should match
        total_processed = sum(len(batch_df) for batch_df in batches)
        assert total_processed == len(variants)
    
    @given(
        miss_fmt=st.sampled_from(["NN", "N", "--", "."]),
        gt_sep=st.sampled_from(["", "/", "|"])
    )
    @settings(max_examples=20)
    def test_format_parameter_consistency(self, miss_fmt, gt_sep):
        """Test that format parameters are applied consistently."""
        converter = GenotypeConverter(miss_fmt=miss_fmt, gt_sep=gt_sep)
        
        # Test missing genotype
        result = converter.convert_genotype("./.", "A", ["T"])
        assert result == miss_fmt
        
        # Test homozygous genotype with single bases
        result = converter.convert_genotype("0/0", "A", ["T"])
        expected = f"A{gt_sep}A" if gt_sep else "AA"
        assert result == expected
        
        # Test heterozygous genotype with single bases
        result = converter.convert_genotype("0/1", "A", ["T"])
        expected = f"A{gt_sep}T" if gt_sep else "AT"
        assert result == expected


class GenotypeConverterStateMachine(RuleBasedStateMachine):
    """Stateful testing for GenotypeConverter operations."""
    
    def __init__(self):
        super().__init__()
        self.converter = GenotypeConverter()
        self.processed_genotypes = []
    
    @rule(
        genotype=valid_genotype_strategy(),
        ref=valid_allele_strategy(),
        alt_alleles=st.lists(valid_allele_strategy(), min_size=1, max_size=2)
    )
    def convert_genotype(self, genotype: str, ref: str, alt_alleles):
        """Convert a genotype and record it."""
        result = self.converter.convert_genotype(genotype, ref, alt_alleles)
        self.processed_genotypes.append((genotype, ref, alt_alleles, result))
    
    @rule()
    def check_stats_consistency(self):
        """Verify that statistics are consistent with processed genotypes."""
        stats = self.converter.get_stats()
        total_processed = len(self.processed_genotypes)
        total_in_stats = stats['total_genotypes']
        
        # Total in stats should match processed genotypes
        assert total_in_stats == total_processed
        
        # Check specific counts
        missing_count = sum(1 for gt, _, _, _ in self.processed_genotypes if gt in ["./.", ".|."])
        assert stats['missing_genotypes'] == missing_count
    
    @rule()
    def get_summary(self):
        """Test summary generation."""
        summary = self.converter.get_summary()
        assert isinstance(summary, str)
        
        if self.processed_genotypes:
            assert "Total:" in summary
        else:
            assert "No genotypes processed" in summary
    
    @rule()
    def reset_stats(self):
        """Reset statistics and verify clean state."""
        self.converter.reset_stats()
        stats = self.converter.get_stats()
        
        # All stats should be zero
        assert all(count == 0 for count in stats.values())
        
        # Clear our tracking too
        self.processed_genotypes.clear()
    
    @invariant()
    def stats_non_negative(self):
        """Statistics should never be negative."""
        stats = self.converter.get_stats()
        assert all(count >= 0 for count in stats.values())
    
    @invariant()
    def results_not_empty(self):
        """Conversion results should never be empty."""
        for _, _, _, result in self.processed_genotypes:
            assert result
            assert isinstance(result, str)


# Run the stateful test
TestGenotypeConverterState = GenotypeConverterStateMachine.TestCase