"""Tests for dual genotype format support (VCF codes vs sequence bases)."""

import tempfile
from pathlib import Path

import pandas as pd
import pytest

from chip.vcf_processor.config import ProcessingConfig, VariantInfo
from chip.vcf_processor.genotype_converter import GenotypeConverter
from chip.vcf_processor.vcf_processor import VCFProcessor


class TestDualGenotypeFormats:
    """Test VCF processor's dual genotype format support."""
    
    def test_vcf_codes_vs_sequence_bases_missing_formats(self):
        """Test that VCF codes use ./. and sequence bases use NN for missing genotypes."""
        # Test VCF codes converter (keep_original_gt=True)
        vcf_converter = GenotypeConverter(
            miss_fmt="./.",
            keep_original_gt=True
        )
        
        # Test sequence bases converter (keep_original_gt=False)  
        seq_converter = GenotypeConverter(
            miss_fmt="NN",
            keep_original_gt=False
        )
        
        # Test missing genotype conversion
        vcf_result = vcf_converter.convert_genotype("./.", "A", ["T"])
        seq_result = seq_converter.convert_genotype("./.", "A", ["T"])
        
        assert vcf_result == "./.", f"VCF codes should use ./. for missing, got {vcf_result}"
        assert seq_result == "NN", f"Sequence bases should use NN for missing, got {seq_result}"
        
        # Test normal genotype conversion
        vcf_normal = vcf_converter.convert_genotype("0/1", "A", ["T"])
        seq_normal = seq_converter.convert_genotype("0/1", "A", ["T"])
        
        assert vcf_normal == "0/1", f"VCF codes should keep original format, got {vcf_normal}"
        assert seq_normal == "AT", f"Sequence bases should convert to bases, got {seq_normal}"
    
    def test_vcf_processor_dual_output_formats(self):
        """Test that VCF processor generates both formats correctly."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test VCF file
            vcf_file = temp_path / "test.vcf"
            vcf_content = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\tsample3
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0\t0/1\t./.
chr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t1/1\t./.\t0/1
"""
            vcf_file.write_text(vcf_content)
            
            # Create target ID file
            target_id_file = temp_path / "targets.txt"
            target_id_file.write_text("chr1_100_A_T\nchr1_200_G_C\n")
            
            # Create output path
            output_file = temp_path / "test_output"
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                compress_output=False,
                quiet=True
            )
            
            # Process VCF
            processor = VCFProcessor(config)
            result = processor.process()
            
            # Check that processing succeeded
            assert len(result.errors) == 0, f"Processing failed with errors: {result.errors}"
            assert result.processed_variants > 0, "Should have processed some variants"
            
            # Check output files exist
            codes_file = output_file.with_suffix(".genotype_codes.tsv")
            bases_file = output_file.with_suffix(".genotype_bases.tsv")
            
            assert codes_file.exists(), f"VCF codes file should exist: {codes_file}"
            assert bases_file.exists(), f"Sequence bases file should exist: {bases_file}"
            
            # Read and verify content
            codes_df = pd.read_csv(codes_file, sep='\t')
            bases_df = pd.read_csv(bases_file, sep='\t')
            
            # Both should have same structure
            assert len(codes_df) == len(bases_df), "Both files should have same number of rows"
            assert list(codes_df.columns) == list(bases_df.columns), "Both files should have same columns"
            
            # Check specific genotype conversions
            # First variant: chr1:100 A>T with genotypes 0/0, 0/1, ./.
            row1_codes = codes_df.iloc[0]
            row1_bases = bases_df.iloc[0]
            
            assert row1_codes['sample1'] == "0/0", "VCF codes should keep original format"
            assert row1_codes['sample2'] == "0/1", "VCF codes should keep original format"
            assert row1_codes['sample3'] == "./.", "VCF codes should use ./. for missing"
            
            assert row1_bases['sample1'] == "AA", "Sequence bases should convert to bases"
            assert row1_bases['sample2'] == "AT", "Sequence bases should convert to bases"
            assert row1_bases['sample3'] == "NN", "Sequence bases should use NN for missing"
            
            # Second variant: chr1:200 G>C with genotypes 1/1, ./., 0/1
            row2_codes = codes_df.iloc[1]
            row2_bases = bases_df.iloc[1]
            
            assert row2_codes['sample1'] == "1/1", "VCF codes should keep original format"
            assert row2_codes['sample2'] == "./.", "VCF codes should use ./. for missing"
            assert row2_codes['sample3'] == "0/1", "VCF codes should keep original format"
            
            assert row2_bases['sample1'] == "CC", "Sequence bases should convert to bases"
            assert row2_bases['sample2'] == "NN", "Sequence bases should use NN for missing"
            assert row2_bases['sample3'] == "GC", "Sequence bases should convert to bases"
    
    def test_file_naming_conventions(self):
        """Test that output files use correct naming conventions."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create minimal test files
            vcf_file = temp_path / "test.vcf"
            vcf_content = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0
"""
            vcf_file.write_text(vcf_content)
            
            target_id_file = temp_path / "targets.txt"
            target_id_file.write_text("chr1_100_A_T\n")
            
            # Test different output base names
            test_cases = [
                "simple_output",
                "output-with-hyphens",
                "output.with.dots",
                "output_with_underscores"
            ]
            
            for base_name in test_cases:
                output_file = temp_path / base_name
                
                config = ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_id_file,
                    output_file=output_file,
                    compress_output=False,
                    quiet=True
                )
                
                processor = VCFProcessor(config)
                result = processor.process()
                
                # Check file naming
                expected_codes = temp_path / f"{base_name}.genotype_codes.tsv"
                expected_bases = temp_path / f"{base_name}.genotype_bases.tsv"
                
                assert expected_codes.exists(), f"Codes file should exist: {expected_codes}"
                assert expected_bases.exists(), f"Bases file should exist: {expected_bases}"
                
                # Cleanup for next iteration
                if expected_codes.exists():
                    expected_codes.unlink()
                if expected_bases.exists():
                    expected_bases.unlink()
    
    def test_compressed_output_naming(self):
        """Test file naming with compression enabled."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create minimal test files
            vcf_file = temp_path / "test.vcf"
            vcf_content = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0
"""
            vcf_file.write_text(vcf_content)
            
            target_id_file = temp_path / "targets.txt"
            target_id_file.write_text("chr1_100_A_T\n")
            
            output_file = temp_path / "compressed_test"
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                compress_output=True,
                quiet=True
            )
            
            processor = VCFProcessor(config)
            result = processor.process()
            
            # Check compressed file naming
            expected_codes = temp_path / "compressed_test.genotype_codes.tsv.gz"
            expected_bases = temp_path / "compressed_test.genotype_bases.tsv.gz"
            
            assert expected_codes.exists(), f"Compressed codes file should exist: {expected_codes}"
            assert expected_bases.exists(), f"Compressed bases file should exist: {expected_bases}"
            
            # Verify files are actually compressed (have .gz extension)
            assert expected_codes.suffix == ".gz", "Codes file should be gzipped"
            assert expected_bases.suffix == ".gz", "Bases file should be gzipped"
    
    def test_multiallelic_variants(self):
        """Test dual format handling of multiallelic variants."""
        # Test converters with multiallelic variants
        vcf_converter = GenotypeConverter(miss_fmt="./.", keep_original_gt=True)
        seq_converter = GenotypeConverter(miss_fmt="NN", keep_original_gt=False)
        
        # Test various multiallelic genotypes
        test_cases = [
            ("0/0", "A", ["T", "G"], "0/0", "AA"),
            ("0/1", "A", ["T", "G"], "0/1", "AT"),
            ("0/2", "A", ["T", "G"], "0/2", "AG"),
            ("1/2", "A", ["T", "G"], "1/2", "TG"),
            ("2/2", "A", ["T", "G"], "2/2", "GG"),
            ("./.", "A", ["T", "G"], "./.", "NN"),
        ]
        
        for genotype, ref, alt, expected_vcf, expected_seq in test_cases:
            vcf_result = vcf_converter.convert_genotype(genotype, ref, alt)
            seq_result = seq_converter.convert_genotype(genotype, ref, alt)
            
            assert vcf_result == expected_vcf, f"VCF format: {genotype} -> expected {expected_vcf}, got {vcf_result}"
            assert seq_result == expected_seq, f"Seq format: {genotype} -> expected {expected_seq}, got {seq_result}"
    
    def test_indel_variants(self):
        """Test dual format handling of insertion/deletion variants."""
        vcf_converter = GenotypeConverter(miss_fmt="./.", keep_original_gt=True)
        seq_converter = GenotypeConverter(miss_fmt="NN", keep_original_gt=False)
        
        # Test indel variants
        test_cases = [
            # Insertion
            ("0/0", "A", ["ATG"], "0/0", "AA"),  # Single base homozygous ref -> AA
            ("0/1", "A", ["ATG"], "0/1", "A/ATG"),  # Heterozygous indels use "/" separator
            ("1/1", "A", ["ATG"], "1/1", "ATG"),  # Multi-base homozygous alt returns single allele
            
            # Deletion
            ("0/0", "ATG", ["A"], "0/0", "ATG"),  # Multi-base homozygous ref returns single allele
            ("0/1", "ATG", ["A"], "0/1", "ATG/A"),
            ("1/1", "ATG", ["A"], "1/1", "AA"),  # Single base homozygous alt -> AA
            
            # Missing
            ("./.", "ATG", ["A"], "./.", "NN"),
        ]
        
        for genotype, ref, alt, expected_vcf, expected_seq in test_cases:
            vcf_result = vcf_converter.convert_genotype(genotype, ref, alt)
            seq_result = seq_converter.convert_genotype(genotype, ref, alt)
            
            assert vcf_result == expected_vcf, f"VCF indel: {genotype} -> expected {expected_vcf}, got {vcf_result}"
            assert seq_result == expected_seq, f"Seq indel: {genotype} -> expected {expected_seq}, got {seq_result}"


class TestGenotypeFormatConsistency:
    """Test consistency between the two genotype formats."""
    
    def test_variant_count_consistency(self):
        """Test that both output files have the same number of variants."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test VCF with multiple variants
            vcf_file = temp_path / "test.vcf"
            vcf_content = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0\t0/1
chr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t1/1\t./.
chr2\t300\t.\tT\tA,G\t60\tPASS\t.\tGT\t0/1\t1/2
"""
            vcf_file.write_text(vcf_content)
            
            target_id_file = temp_path / "targets.txt"
            target_id_file.write_text("chr1_100_A_T\nchr1_200_G_C\nchr2_300_T_A\n")
            
            output_file = temp_path / "consistency_test"
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                compress_output=False,
                quiet=True
            )
            
            processor = VCFProcessor(config)
            result = processor.process()
            
            # Read both output files
            codes_file = output_file.with_suffix(".genotype_codes.tsv")
            bases_file = output_file.with_suffix(".genotype_bases.tsv")
            
            codes_df = pd.read_csv(codes_file, sep='\t')
            bases_df = pd.read_csv(bases_file, sep='\t')
            
            # Should have same dimensions
            assert codes_df.shape == bases_df.shape, "Both files should have same dimensions"
            
            # Should have same variant positions
            assert list(codes_df['CHROM']) == list(bases_df['CHROM']), "CHROM columns should match"
            assert list(codes_df['POS']) == list(bases_df['POS']), "POS columns should match"
            assert list(codes_df['REF']) == list(bases_df['REF']), "REF columns should match"
            assert list(codes_df['ALT']) == list(bases_df['ALT']), "ALT columns should match"
    
    def test_sample_consistency(self):
        """Test that both files have the same sample columns."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test VCF with many samples
            vcf_file = temp_path / "test.vcf"
            samples = ["sample1", "sample2", "sample3", "sample4", "sample5"]
            header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)
            genotypes = "\t".join(["0/0", "0/1", "1/1", "./.", "0/1"])
            
            vcf_content = f"""##fileformat=VCFv4.2
{header}
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t{genotypes}
"""
            vcf_file.write_text(vcf_content)
            
            target_id_file = temp_path / "targets.txt"
            target_id_file.write_text("chr1_100_A_T\n")
            
            output_file = temp_path / "sample_test"
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                compress_output=False,
                quiet=True
            )
            
            processor = VCFProcessor(config)
            result = processor.process()
            
            # Read both output files
            codes_file = output_file.with_suffix(".genotype_codes.tsv")
            bases_file = output_file.with_suffix(".genotype_bases.tsv")
            
            codes_df = pd.read_csv(codes_file, sep='\t')
            bases_df = pd.read_csv(bases_file, sep='\t')
            
            # Should have same column names
            assert list(codes_df.columns) == list(bases_df.columns), "Column names should match"
            
            # Sample columns should be present
            for sample in samples:
                assert sample in codes_df.columns, f"Sample {sample} should be in codes file"
                assert sample in bases_df.columns, f"Sample {sample} should be in bases file"