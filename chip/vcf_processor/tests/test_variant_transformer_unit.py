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