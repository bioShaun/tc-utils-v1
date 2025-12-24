"""Unit tests for genotype conversion functionality."""

import pandas as pd
import pytest
from ..genotype_converter import GenotypeConverter, BatchProcessor
from ..config import VariantInfo


class TestGenotypeConverterUnit:
    """Unit tests for specific genotype conversion cases."""
    
    def test_missing_genotype_conversion(self):
        """Test conversion of missing genotypes."""
        converter = GenotypeConverter(miss_fmt="NN")
        
        # Test various missing formats
        assert converter.convert_genotype("./.", "A", ["T"]) == "NN"
        assert converter.convert_genotype(".|.", "A", ["T"]) == "NN"
        
        # Test custom missing format
        converter_custom = GenotypeConverter(miss_fmt="--")
        assert converter_custom.convert_genotype("./.", "A", ["T"]) == "--"
        
        # Verify stats
        stats = converter.get_stats()
        assert stats['missing_genotypes'] == 2
        assert stats['total_genotypes'] == 2
    
    def test_homozygous_ref_conversion(self):
        """Test conversion of homozygous reference genotypes."""
        converter = GenotypeConverter(gt_sep="")
        
        # Single base reference
        assert converter.convert_genotype("0/0", "A", ["T"]) == "AA"
        assert converter.convert_genotype("0|0", "G", ["C"]) == "GG"
        
        # Multi-base reference
        assert converter.convert_genotype("0/0", "ATG", ["C"]) == "ATG"
        
        # Verify stats
        stats = converter.get_stats()
        assert stats['homozygous_ref'] == 3
        assert stats['total_genotypes'] == 3
    
    def test_homozygous_alt_conversion(self):
        """Test conversion of homozygous alternative genotypes."""
        converter = GenotypeConverter(gt_sep="")
        
        # Single base alternative
        assert converter.convert_genotype("1/1", "A", ["T"]) == "TT"
        assert converter.convert_genotype("2/2", "A", ["T", "G"]) == "GG"
        
        # Multi-base alternative
        assert converter.convert_genotype("1/1", "A", ["ATG"]) == "ATG"
        
        # Verify stats
        stats = converter.get_stats()
        assert stats['homozygous_alt'] == 3
        assert stats['total_genotypes'] == 3
    
    def test_heterozygous_conversion(self):
        """Test conversion of heterozygous genotypes."""
        converter = GenotypeConverter(gt_sep="")
        
        # Single base heterozygous
        assert converter.convert_genotype("0/1", "A", ["T"]) == "AT"
        assert converter.convert_genotype("1/0", "A", ["T"]) == "TA"
        assert converter.convert_genotype("0|1", "G", ["C"]) == "GC"
        
        # Multi-base heterozygous (should use "/" separator)
        assert converter.convert_genotype("0/1", "ATG", ["C"]) == "ATG/C"
        assert converter.convert_genotype("1/2", "A", ["ATG", "C"]) == "ATG/C"
        
        # Verify stats
        stats = converter.get_stats()
        assert stats['heterozygous'] == 5
        assert stats['total_genotypes'] == 5
    
    def test_custom_separator(self):
        """Test custom genotype separator."""
        converter = GenotypeConverter(gt_sep="|")
        
        # Single base with custom separator
        assert converter.convert_genotype("0/0", "A", ["T"]) == "A|A"
        assert converter.convert_genotype("0/1", "A", ["T"]) == "A|T"
        
        # Multi-base should still use "/" for heterozygous
        assert converter.convert_genotype("0/1", "ATG", ["C"]) == "ATG/C"
    
    def test_error_handling(self):
        """Test error handling for invalid genotypes."""
        converter = GenotypeConverter(miss_fmt="NN")
        
        # Invalid allele index
        result = converter.convert_genotype("5/5", "A", ["T"])
        assert result == "NN"
        
        # Invalid format
        result = converter.convert_genotype("invalid", "A", ["T"])
        assert result == "NN"
        
        # Verify error stats
        stats = converter.get_stats()
        assert stats['errors'] == 2
    
    def test_variant_genotype_conversion(self):
        """Test conversion of all genotypes for a variant."""
        converter = GenotypeConverter()
        
        variant = VariantInfo(
            chrom="chr1",
            pos=100,
            ref="A",
            alt=["T", "G"],
            variant_id="chr1_100",
            genotypes=["0/0", "0/1", "1/1", "0/2", "./.", "1/2"]
        )
        
        sample_names = ["sample1", "sample2", "sample3", "sample4", "sample5", "sample6"]
        
        result = converter.convert_variant_genotypes(variant, sample_names)
        
        assert result["sample1"] == "AA"  # 0/0
        assert result["sample2"] == "AT"  # 0/1
        assert result["sample3"] == "TT"  # 1/1
        assert result["sample4"] == "AG"  # 0/2
        assert result["sample5"] == "NN"  # ./.
        assert result["sample6"] == "TG"  # 1/2
    
    def test_batch_conversion(self):
        """Test batch conversion to DataFrame."""
        converter = GenotypeConverter()
        
        variants = [
            VariantInfo(
                chrom="chr1", pos=100, ref="A", alt=["T"],
                variant_id="chr1_100", genotypes=["0/0", "0/1"]
            ),
            VariantInfo(
                chrom="chr1", pos=200, ref="G", alt=["C"],
                variant_id="chr1_200", genotypes=["1/1", "./."]
            )
        ]
        
        sample_names = ["sample1", "sample2"]
        
        df = converter.convert_batch(variants, sample_names)
        
        # Check DataFrame structure
        assert len(df) == 2
        assert list(df.columns) == ["CHROM", "POS", "REF", "ALT", "sample1", "sample2"]
        
        # Check first variant
        assert df.iloc[0]["CHROM"] == "chr1"
        assert df.iloc[0]["POS"] == 100
        assert df.iloc[0]["REF"] == "A"
        assert df.iloc[0]["ALT"] == "T"
        assert df.iloc[0]["sample1"] == "AA"
        assert df.iloc[0]["sample2"] == "AT"
        
        # Check second variant
        assert df.iloc[1]["sample1"] == "CC"
        assert df.iloc[1]["sample2"] == "NN"
    
    def test_dataframe_processing(self):
        """Test processing of existing DataFrame."""
        converter = GenotypeConverter()
        
        # Create test DataFrame
        df = pd.DataFrame({
            'CHROM': ['chr1', 'chr2'],
            'POS': [100, 200],
            'REF': ['A', 'G'],
            'ALT': ['T', 'C'],
            'sample1': ['0/0', '1/1'],
            'sample2': ['0/1', './.']
        })
        
        result_df = converter.process_dataframe(df)
        
        # Check conversions
        assert result_df.loc[0, 'sample1'] == "AA"
        assert result_df.loc[0, 'sample2'] == "AT"
        assert result_df.loc[1, 'sample1'] == "CC"
        assert result_df.loc[1, 'sample2'] == "NN"
        
        # Original DataFrame should be unchanged
        assert df.loc[0, 'sample1'] == "0/0"
    
    def test_statistics_and_summary(self):
        """Test statistics tracking and summary generation."""
        converter = GenotypeConverter()
        
        # Process various genotype types
        converter.convert_genotype("0/0", "A", ["T"])  # Homozygous ref
        converter.convert_genotype("1/1", "A", ["T"])  # Homozygous alt
        converter.convert_genotype("0/1", "A", ["T"])  # Heterozygous
        converter.convert_genotype("./.", "A", ["T"])  # Missing
        converter.convert_genotype("invalid", "A", ["T"])  # Error
        
        # Check stats
        stats = converter.get_stats()
        assert stats['total_genotypes'] == 5
        assert stats['homozygous_ref'] == 1
        assert stats['homozygous_alt'] == 1
        assert stats['heterozygous'] == 1
        assert stats['missing_genotypes'] == 1
        assert stats['errors'] == 1
        
        # Check summary
        summary = converter.get_summary()
        assert "Total: 5" in summary
        assert "Homozygous Ref: 1 (20.0%)" in summary
        assert "Missing: 1 (20.0%)" in summary
        
        # Reset and check
        converter.reset_stats()
        stats = converter.get_stats()
        assert all(count == 0 for count in stats.values())
        
        summary = converter.get_summary()
        assert "No genotypes processed" in summary


class TestBatchProcessor:
    """Unit tests for BatchProcessor class."""
    
    def test_batch_processing(self):
        """Test basic batch processing functionality."""
        processor = BatchProcessor(batch_size=2)
        converter = GenotypeConverter()
        
        variants = [
            VariantInfo(chrom="chr1", pos=i, ref="A", alt=["T"], 
                       variant_id=f"chr1_{i}", genotypes=["0/0", "0/1"])
            for i in range(100, 106)  # 6 variants
        ]
        
        sample_names = ["sample1", "sample2"]
        
        batches = list(processor.process_variants(iter(variants), sample_names, converter))
        
        # Should have 3 batches (2+2+2)
        assert len(batches) == 3
        assert len(batches[0]) == 2
        assert len(batches[1]) == 2
        assert len(batches[2]) == 2
        
        # Check processed count
        assert processor.get_processed_count() == 6
    
    def test_partial_last_batch(self):
        """Test handling of partial last batch."""
        processor = BatchProcessor(batch_size=3)
        converter = GenotypeConverter()
        
        variants = [
            VariantInfo(chrom="chr1", pos=i, ref="A", alt=["T"], 
                       variant_id=f"chr1_{i}", genotypes=["0/0"])
            for i in range(100, 105)  # 5 variants
        ]
        
        sample_names = ["sample1"]
        
        batches = list(processor.process_variants(iter(variants), sample_names, converter))
        
        # Should have 2 batches (3+2)
        assert len(batches) == 2
        assert len(batches[0]) == 3
        assert len(batches[1]) == 2
        
        assert processor.get_processed_count() == 5
    
    def test_progress_callback(self):
        """Test progress callback functionality."""
        progress_calls = []
        
        def progress_callback(total_processed, batch_size):
            progress_calls.append((total_processed, batch_size))
        
        processor = BatchProcessor(batch_size=2, progress_callback=progress_callback)
        converter = GenotypeConverter()
        
        variants = [
            VariantInfo(chrom="chr1", pos=i, ref="A", alt=["T"], 
                       variant_id=f"chr1_{i}", genotypes=["0/0"])
            for i in range(100, 104)  # 4 variants
        ]
        
        sample_names = ["sample1"]
        
        list(processor.process_variants(iter(variants), sample_names, converter))
        
        # Should have 2 progress calls
        assert len(progress_calls) == 2
        assert progress_calls[0] == (2, 2)  # First batch
        assert progress_calls[1] == (4, 2)  # Second batch
    
    def test_reset_count(self):
        """Test resetting processed count."""
        processor = BatchProcessor(batch_size=1)
        converter = GenotypeConverter()
        
        variants = [
            VariantInfo(chrom="chr1", pos=100, ref="A", alt=["T"], 
                       variant_id="chr1_100", genotypes=["0/0"])
        ]
        
        sample_names = ["sample1"]
        
        list(processor.process_variants(iter(variants), sample_names, converter))
        assert processor.get_processed_count() == 1
        
        processor.reset_count()
        assert processor.get_processed_count() == 0