"""Unit tests for variant transformation edge cases."""

import pytest
from ..variant_transformer import VariantTransformer, validate_variant_alleles
from ..config import VariantInfo


class TestVariantTransformerUnit:
    """Unit tests for specific variant transformation cases."""
    
    def test_snp_transformation(self):
        """Test SNP transformation cases."""
        transformer = VariantTransformer()
        
        # Simple SNPs
        assert transformer.transform_one_alt("A", "T") == "T"
        assert transformer.transform_one_alt("G", "C") == "C"
        assert transformer.transform_one_alt("T", "A") == "A"
        assert transformer.transform_one_alt("C", "G") == "G"
        
        # Verify stats
        stats = transformer.get_stats()
        assert stats['snps'] == 4
    
    def test_insertion_transformation(self):
        """Test insertion transformation cases."""
        transformer = VariantTransformer()
        
        # Simple insertions
        assert transformer.transform_one_alt("A", "AT") == "insT"
        assert transformer.transform_one_alt("G", "GCA") == "insCA"
        assert transformer.transform_one_alt("T", "TAAA") == "insAAA"
        
        # Complex insertions with common prefix/suffix
        assert transformer.transform_one_alt("ATG", "ATCG") == "insC"
        
        # Verify stats
        stats = transformer.get_stats()
        assert stats['insertions'] == 4
    
    def test_deletion_transformation(self):
        """Test deletion transformation cases."""
        transformer = VariantTransformer()
        
        # Simple deletions
        assert transformer.transform_one_alt("AT", "A") == "delT"
        assert transformer.transform_one_alt("GCA", "G") == "delCA"
        assert transformer.transform_one_alt("TAAA", "T") == "delAAA"
        
        # Symbolic deletion
        assert transformer.transform_one_alt("ATG", "*") == "delATG"
        
        # Complex deletions with common prefix/suffix
        assert transformer.transform_one_alt("ATCG", "ATG") == "delC"
        
        # Verify stats
        stats = transformer.get_stats()
        assert stats['deletions'] == 5
    
    def test_mnp_transformation(self):
        """Test multi-nucleotide polymorphism transformation."""
        transformer = VariantTransformer()
        
        # MNPs (same length, different sequence)
        assert transformer.transform_one_alt("AT", "GC") == "GC"
        assert transformer.transform_one_alt("GCA", "TTT") == "TTT"
        assert transformer.transform_one_alt("AAAA", "TTTT") == "TTTT"
        
        # Verify stats
        stats = transformer.get_stats()
        assert stats['mnps'] == 3
    
    def test_complex_variant_strict_mode(self):
        """Test complex variant handling in strict mode."""
        transformer = VariantTransformer(strict_mode=True)
        
        # Complex variants should raise exceptions
        with pytest.raises(ValueError, match="Complex variant"):
            transformer.transform_one_alt("AT", "GCA")  # Different lengths, not simple indel
        
        with pytest.raises(ValueError, match="Complex variant"):
            transformer.transform_one_alt("ATCG", "GG")  # Complex delins
    
    def test_complex_variant_lenient_mode(self):
        """Test complex variant handling in lenient mode."""
        transformer = VariantTransformer(strict_mode=False)
        
        # Complex variants should return original alt
        assert transformer.transform_one_alt("AT", "GCA") == "GCA"
        assert transformer.transform_one_alt("ATCG", "GG") == "GG"
        
        # Verify stats
        stats = transformer.get_stats()
        assert stats['complex'] == 2
    
    def test_error_handling(self):
        """Test error handling for invalid inputs."""
        transformer_strict = VariantTransformer(strict_mode=True)
        transformer_lenient = VariantTransformer(strict_mode=False)
        
        # Empty inputs should raise in strict mode
        with pytest.raises(ValueError, match="empty"):
            transformer_strict.transform_one_alt("", "A")
        
        with pytest.raises(ValueError, match="empty"):
            transformer_strict.transform_one_alt("A", "")
        
        # Lenient mode should handle gracefully
        assert transformer_lenient.transform_one_alt("", "A") == "A"
        assert transformer_lenient.transform_one_alt("A", "") == ""
        
        # Verify error stats
        stats = transformer_lenient.get_stats()
        assert stats['errors'] == 2
    
    def test_prefix_suffix_removal(self):
        """Test proper handling of common prefix and suffix removal."""
        transformer = VariantTransformer()
        
        # Common prefix only
        assert transformer.transform_one_alt("ATCG", "ATTT") == "TT"  # MNP after prefix removal
        
        # Common suffix only - ATCG vs TTCG has common suffix TCG, leaving A vs T (SNP)
        assert transformer.transform_one_alt("ATCG", "TTCG") == "T"  # SNP after suffix removal
        
        # Both common prefix and suffix
        assert transformer.transform_one_alt("ATCGAA", "ATTTAA") == "TT"  # MNP in middle
        
        # Insertion with common prefix/suffix
        assert transformer.transform_one_alt("ATCG", "ATCCG") == "insC"
        
        # Deletion with common prefix/suffix
        assert transformer.transform_one_alt("ATCCG", "ATCG") == "delC"
    
    def test_multiple_alt_alleles(self):
        """Test handling of multiple alternative alleles."""
        transformer = VariantTransformer(strict_mode=False)  # Use lenient mode
        
        # Multiple alts
        alt_list = ["T", "G", "C"]
        results = transformer.transform_alt_alleles("A", alt_list)
        assert results == ["T", "G", "C"]  # All SNPs
        
        # Mixed variant types
        alt_list = ["T", "AT", ""]  # SNP, insertion, empty (error case)
        results = transformer.transform_alt_alleles("A", alt_list)
        assert results[0] == "T"  # SNP
        assert results[1] == "insT"  # Insertion
        assert results[2] == ""  # Empty (error case in lenient mode)
        
        # String format
        result_string = transformer.transform_alt_string("A", "T,AT")
        assert result_string == "T,insT"
    
    def test_variant_info_transformation(self):
        """Test transformation of VariantInfo objects."""
        transformer = VariantTransformer(strict_mode=False)  # Use lenient mode
        
        # Create test variant with simpler cases
        variant = VariantInfo(
            chrom="chr1",
            pos=100,
            ref="A",  # Single base ref
            alt=["T", "AT"],  # SNP and insertion
            variant_id="chr1_100",
            genotypes=["0/1", "1/2"]
        )
        
        # Transform
        transformed = transformer.transform_variant(variant)
        
        # Check results
        assert transformed.chrom == "chr1"
        assert transformed.pos == 100
        assert transformed.ref == "A"  # Single base
        assert transformed.alt == ["T", "insT"]  # Transformed alts
        assert transformed.variant_id == "chr1_100"
        assert transformed.genotypes == ["0/1", "1/2"]
    
    def test_statistics_and_summary(self):
        """Test statistics tracking and summary generation."""
        transformer = VariantTransformer()
        
        # Process various variant types
        transformer.transform_one_alt("A", "T")  # SNP
        transformer.transform_one_alt("A", "AT")  # Insertion
        transformer.transform_one_alt("AT", "A")  # Deletion
        transformer.transform_one_alt("AT", "GC")  # MNP
        
        # Check stats
        stats = transformer.get_stats()
        assert stats['snps'] == 1
        assert stats['insertions'] == 1
        assert stats['deletions'] == 1
        assert stats['mnps'] == 1
        assert stats['complex'] == 0
        assert stats['errors'] == 0
        
        # Check summary
        summary = transformer.get_summary()
        assert "Total: 4" in summary
        assert "SNPs: 1 (25.0%)" in summary
        assert "Insertions: 1 (25.0%)" in summary
        
        # Reset and check
        transformer.reset_stats()
        stats = transformer.get_stats()
        assert all(count == 0 for count in stats.values())
        
        summary = transformer.get_summary()
        assert "No variants processed" in summary
    
    def test_validation_function(self):
        """Test the validation helper function thoroughly."""
        # Valid cases
        assert validate_variant_alleles("A", ["T"]) == []
        assert validate_variant_alleles("ATG", ["A", "T", "G"]) == []
        assert validate_variant_alleles("C", ["*"]) == []  # Symbolic
        assert validate_variant_alleles("ATCGN", ["ATCGN"]) == []  # With N
        
        # Invalid reference
        issues = validate_variant_alleles("", ["T"])
        assert len(issues) == 1
        assert "empty" in issues[0].lower()
        
        issues = validate_variant_alleles("AX", ["T"])
        assert len(issues) == 1
        assert "invalid characters" in issues[0].lower()
        
        # Invalid alternatives
        issues = validate_variant_alleles("A", [])
        assert len(issues) == 1
        assert "no alternative" in issues[0].lower()
        
        issues = validate_variant_alleles("A", [""])
        assert len(issues) == 1
        assert "empty" in issues[0].lower()
        
        issues = validate_variant_alleles("A", ["X"])
        assert len(issues) == 1
        assert "invalid characters" in issues[0].lower()
        
        # Multiple issues
        issues = validate_variant_alleles("", ["", "X"])
        assert len(issues) == 3  # Empty ref, empty alt1, invalid alt2


class TestVariantTypeClassification:
    """Unit tests for variant type classification edge cases."""
    
    def test_classify_variant_type_snp_cases(self):
        """Test SNP classification edge cases."""
        # Standard SNPs
        assert VariantTransformer.classify_variant_type("A", "T") == "SNP"
        assert VariantTransformer.classify_variant_type("G", "C") == "SNP"
        assert VariantTransformer.classify_variant_type("T", "A") == "SNP"
        assert VariantTransformer.classify_variant_type("C", "G") == "SNP"
        
        # Case sensitivity (should work with any case)
        assert VariantTransformer.classify_variant_type("a", "t") == "SNP"
        assert VariantTransformer.classify_variant_type("G", "c") == "SNP"
    
    def test_classify_variant_type_indel_cases(self):
        """Test INDEL classification edge cases."""
        # Simple insertions
        assert VariantTransformer.classify_variant_type("A", "AT") == "INDEL"
        assert VariantTransformer.classify_variant_type("G", "GCA") == "INDEL"
        assert VariantTransformer.classify_variant_type("T", "TAAAAAA") == "INDEL"
        
        # Simple deletions
        assert VariantTransformer.classify_variant_type("AT", "A") == "INDEL"
        assert VariantTransformer.classify_variant_type("GCA", "G") == "INDEL"
        assert VariantTransformer.classify_variant_type("TAAAAAA", "T") == "INDEL"
        
        # Symbolic deletion
        assert VariantTransformer.classify_variant_type("A", "*") == "INDEL"
        assert VariantTransformer.classify_variant_type("ATG", "*") == "INDEL"
        
        # Complex indels (very different lengths)
        assert VariantTransformer.classify_variant_type("A", "ATCGATCGATCG") == "INDEL"
        assert VariantTransformer.classify_variant_type("ATCGATCGATCG", "A") == "INDEL"
    
    def test_classify_variant_type_mnp_cases(self):
        """Test MNP classification edge cases."""
        # Standard MNPs
        assert VariantTransformer.classify_variant_type("AT", "GC") == "MNP"
        assert VariantTransformer.classify_variant_type("GCA", "TTT") == "MNP"
        assert VariantTransformer.classify_variant_type("AAAA", "TTTT") == "MNP"
        
        # Long MNPs
        assert VariantTransformer.classify_variant_type("ATCGATCG", "GCTAGCTA") == "MNP"
        
        # MNPs with repeated sequences
        assert VariantTransformer.classify_variant_type("AAAA", "GGGG") == "MNP"
        assert VariantTransformer.classify_variant_type("ATATATAT", "GCGCGCGC") == "MNP"
    
    def test_classify_variant_type_ref_cases(self):
        """Test REF classification edge cases."""
        # Identical sequences
        assert VariantTransformer.classify_variant_type("A", "A") == "REF"
        assert VariantTransformer.classify_variant_type("ATG", "ATG") == "REF"
        assert VariantTransformer.classify_variant_type("ATCGATCGATCG", "ATCGATCGATCG") == "REF"
        
        # Case sensitivity
        assert VariantTransformer.classify_variant_type("A", "a") == "REF"
        assert VariantTransformer.classify_variant_type("ATG", "atg") == "REF"
    
    def test_classify_variant_type_edge_cases(self):
        """Test edge cases and error conditions."""
        # Empty inputs
        assert VariantTransformer.classify_variant_type("", "A") == "UNKNOWN"
        assert VariantTransformer.classify_variant_type("A", "") == "UNKNOWN"
        assert VariantTransformer.classify_variant_type("", "") == "UNKNOWN"
        
        # Invalid nucleotides (should still classify by length rules)
        assert VariantTransformer.classify_variant_type("X", "Y") == "SNP"  # Same length, different
        assert VariantTransformer.classify_variant_type("X", "XY") == "INDEL"  # Different length
        assert VariantTransformer.classify_variant_type("XY", "ZW") == "MNP"  # Same length >1
        assert VariantTransformer.classify_variant_type("X", "X") == "REF"  # Identical
    
    def test_classify_multi_allelic_types_basic(self):
        """Test basic multi-allelic type combination."""
        # Single type cases
        assert VariantTransformer.classify_multi_allelic_types("A", ["T"]) == "SNP"
        assert VariantTransformer.classify_multi_allelic_types("A", ["AT"]) == "INDEL"
        assert VariantTransformer.classify_multi_allelic_types("AT", ["GC"]) == "MNP"
        assert VariantTransformer.classify_multi_allelic_types("A", ["A"]) == "REF"
        
        # Multiple types - should be sorted
        assert VariantTransformer.classify_multi_allelic_types("A", ["T", "AT"]) == "SNP|INDEL"
        assert VariantTransformer.classify_multi_allelic_types("A", ["AT", "T"]) == "SNP|INDEL"  # Order shouldn't matter
        
        # All types - need to include a real MNP (same length >1)
        assert VariantTransformer.classify_multi_allelic_types("AT", ["GC", "A", "ATG", "AT"]) == "INDEL|MNP|REF"
        
        # Reverse order input - should still be sorted correctly  
        assert VariantTransformer.classify_multi_allelic_types("AT", ["AT", "ATG", "A", "GC"]) == "INDEL|MNP|REF"
    
    def test_classify_multi_allelic_types_edge_cases(self):
        """Test edge cases for multi-allelic classification."""
        # Empty alt list
        assert VariantTransformer.classify_multi_allelic_types("A", []) == "REF"
        
        # All same type
        assert VariantTransformer.classify_multi_allelic_types("A", ["T", "G", "C"]) == "SNP"
        assert VariantTransformer.classify_multi_allelic_types("A", ["AT", "AG", "AC"]) == "INDEL"
        
        # With symbolic deletion
        assert VariantTransformer.classify_multi_allelic_types("A", ["T", "*"]) == "SNP|INDEL"
        
        # With unknown/invalid types (should be filtered out)
        assert VariantTransformer.classify_multi_allelic_types("", ["T", "G"]) == "UNKNOWN"
        
        # Duplicates should be deduplicated
        assert VariantTransformer.classify_multi_allelic_types("A", ["T", "T", "AT", "AT"]) == "SNP|INDEL"
    
    def test_get_variant_type_column(self):
        """Test variant type column generation."""
        transformer = VariantTransformer()
        
        # Create test variants
        variants = [
            VariantInfo(chrom="chr1", pos=100, ref="A", alt=["T"], variant_id="chr1_100", genotypes=["0/1"]),
            VariantInfo(chrom="chr1", pos=200, ref="A", alt=["AT"], variant_id="chr1_200", genotypes=["0/1"]),
            VariantInfo(chrom="chr1", pos=300, ref="AT", alt=["GC"], variant_id="chr1_300", genotypes=["0/1"]),
            VariantInfo(chrom="chr1", pos=400, ref="A", alt=["A"], variant_id="chr1_400", genotypes=["0/1"]),
            VariantInfo(chrom="chr1", pos=500, ref="A", alt=["T", "AT"], variant_id="chr1_500", genotypes=["0/1"]),
        ]
        
        # Get variant type column
        variant_types = transformer.get_variant_type_column(variants)
        
        # Check results
        assert len(variant_types) == 5
        assert variant_types[0] == "SNP"
        assert variant_types[1] == "INDEL"
        assert variant_types[2] == "MNP"
        assert variant_types[3] == "REF"
        assert variant_types[4] == "SNP|INDEL"
    
    def test_variant_type_with_complex_cases(self):
        """Test variant type classification with complex real-world cases."""
        # Complex structural variants
        assert VariantTransformer.classify_variant_type("ATCGATCG", "A") == "INDEL"  # Large deletion
        assert VariantTransformer.classify_variant_type("A", "ATCGATCGATCGATCG") == "INDEL"  # Large insertion
        
        # Microsatellite-like variants
        assert VariantTransformer.classify_variant_type("ATATAT", "ATATATAT") == "INDEL"  # Repeat expansion
        assert VariantTransformer.classify_variant_type("ATATATAT", "ATATAT") == "INDEL"  # Repeat contraction
        
        # Complex substitutions
        assert VariantTransformer.classify_variant_type("ATCG", "GCTA") == "MNP"  # Complex MNP
        
        # Multi-allelic with all types
        complex_alts = ["T", "ATCG", "GC", "A", "*"]  # SNP, INDEL, MNP, REF, symbolic
        result = VariantTransformer.classify_multi_allelic_types("A", complex_alts)
        assert result == "SNP|INDEL|REF"  # MNP not possible with ref="A"
        
        # Real-world multi-allelic example
        real_alts = ["G", "GA", "GAA"]  # SNP and two different insertions
        result = VariantTransformer.classify_multi_allelic_types("A", real_alts)
        assert result == "SNP|INDEL"