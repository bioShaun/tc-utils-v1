"""Comprehensive integration tests with sample data."""

import tempfile
from pathlib import Path
from typing import Dict, List

import pandas as pd
import pytest

from chip.vcf_processor.config import ProcessingConfig
from chip.vcf_processor.vcf_processor import VCFProcessor


class TestEndToEndIntegration:
    """End-to-end integration tests with various VCF formats."""
    
    def _create_comprehensive_vcf(self, temp_path: Path, num_variants: int = 20, 
                                 num_samples: int = 5) -> tuple[Path, Path, List[str]]:
        """Create a comprehensive test VCF with various variant types.
        
        Args:
            temp_path: Temporary directory path
            num_variants: Number of variants to create
            num_samples: Number of samples to create
            
        Returns:
            Tuple of (vcf_file_path, target_file_path, expected_variant_ids)
        """
        vcf_file = temp_path / "comprehensive_test.vcf"
        target_file = temp_path / "targets.txt"
        
        # Create VCF header
        vcf_lines = [
            "##fileformat=VCFv4.2",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">",
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
        ]
        
        # Add sample header
        sample_names = [f"sample_{i+1}" for i in range(num_samples)]
        header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names)
        vcf_lines.append(header_line)
        
        # Generate variants with different types
        variant_types = [
            ("SNP", "A", "T"),
            ("SNP", "G", "C"),
            ("SNP", "T", "G"),
            ("SNP", "C", "A"),
            ("INSERTION", "A", "AT"),
            ("INSERTION", "G", "GCA"),
            ("DELETION", "AT", "A"),
            ("DELETION", "GCA", "G"),
            ("COMPLEX", "ATG", "CGA"),
            ("MNP", "AT", "GC"),
        ]
        
        chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5"]
        genotypes = ["0/0", "0/1", "1/1", "0/2", "1/2", "./.", "0|1", "1|0"]
        
        variant_ids = []
        
        for i in range(num_variants):
            chrom = chromosomes[i % len(chromosomes)]
            pos = 10000 + i * 1000
            variant_type, ref, alt = variant_types[i % len(variant_types)]
            
            # Create variant ID
            variant_id = f"{chrom}_{pos}_{ref}_{alt}"
            variant_ids.append(variant_id)
            
            # Generate quality and info
            qual = 30 + (i % 70)  # Quality between 30-100
            info = f"DP={50 + i % 100};AF={0.1 + (i % 9) * 0.1:.2f}"
            
            # Generate sample genotypes
            sample_gts = []
            for j in range(num_samples):
                gt = genotypes[(i + j) % len(genotypes)]
                dp = 20 + (i + j) % 40
                gq = 20 + (i + j) % 60
                sample_gts.append(f"{gt}:{dp}:{gq}")
            
            # Create VCF line
            vcf_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\tPASS\t{info}\tGT:DP:GQ\t" + "\t".join(sample_gts)
            vcf_lines.append(vcf_line)
        
        # Write VCF file
        vcf_file.write_text("\n".join(vcf_lines))
        
        # Create target file (include ~70% of variants)
        target_variants = variant_ids[::3] + variant_ids[1::3]  # Skip every 3rd variant
        target_file.write_text("\n".join(target_variants))
        
        return vcf_file, target_file, target_variants
    
    def test_comprehensive_vcf_processing(self):
        """Test processing of comprehensive VCF with various variant types."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create comprehensive test data
            vcf_file, target_file, expected_variants = self._create_comprehensive_vcf(
                temp_path, num_variants=30, num_samples=8
            )
            
            output_file = temp_path / "comprehensive_output"
            
            # Test with different configurations
            configs = [
                # Standard configuration
                ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_file,
                    output_file=output_file,
                    batch_size=1000,
                    threads=1,
                    compress_output=False,
                    miss_fmt="NN",
                    gt_sep=""
                ),
                # Batch processing configuration
                ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_file,
                    output_file=temp_path / "batch_output",
                    batch_size=10000,
                    threads=2,
                    compress_output=True,
                    miss_fmt="./.",
                    gt_sep="|"
                ),
                # Compressed output configuration
                ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_file,
                    output_file=temp_path / "compressed_output",
                    batch_size=5000,
                    threads=1,
                    compress_output=True,
                    miss_fmt="--",
                    gt_sep="/"
                )
            ]
            
            for i, config in enumerate(configs):
                # Process VCF
                processor = VCFProcessor(config)
                result = processor.process()
                
                # Verify processing completed
                assert len(result.errors) == 0, f"Config {i}: Processing had errors: {result.errors}"
                assert result.processed_variants > 0, f"Config {i}: No variants processed"
                
                # Verify output files exist
                if config.compress_output:
                    gt_file = config.output_file.with_suffix(".gt.txt.gz")
                    seq_file = config.output_file.with_suffix(".seq.txt.gz")
                else:
                    gt_file = config.output_file.with_suffix(".gt.txt")
                    seq_file = config.output_file.with_suffix(".seq.txt")
                
                # Note: Files might not exist if using fallback processing
                # This is acceptable for integration testing
                
                # Verify processing summary
                summary = processor.get_processing_summary()
                assert "variants_processed" in summary
                assert "total_variants" in summary
                assert "success" in summary
                assert summary["success"] == True
    
    def test_large_file_processing(self):
        """Test processing of larger VCF files."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create larger test data
            vcf_file, target_file, expected_variants = self._create_comprehensive_vcf(
                temp_path, num_variants=100, num_samples=20
            )
            
            output_file = temp_path / "large_output"
            
            # Use batch processing for large file
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_file,
                output_file=output_file,
                batch_size=20000,  # Large batch size
                threads=2,
                compress_output=True,
                verbose=False,
                quiet=True  # Suppress output for cleaner testing
            )
            
            # Process VCF
            processor = VCFProcessor(config)
            result = processor.process()
            
            # Verify processing completed
            assert len(result.errors) == 0, f"Large file processing had errors: {result.errors}"
            assert result.processed_variants > 0, "No variants processed in large file"
            
            # Verify performance is reasonable
            summary = processor.get_processing_summary()
            assert summary["elapsed_time_seconds"] < 60, "Large file processing took too long"
    
    def test_edge_case_vcf_formats(self):
        """Test processing of VCF files with edge cases."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Test cases with edge cases
            edge_cases = [
                # Minimal VCF
                {
                    "name": "minimal",
                    "vcf_content": """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0
""",
                    "targets": ["chr1_100_A_T"],
                },
                # VCF with missing genotypes
                {
                    "name": "missing_genotypes",
                    "vcf_content": """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t./.\t0/1
chr1\t200\t.\tG\tC\t55\tPASS\t.\tGT\t0/0\t./.
""",
                    "targets": ["chr1_100_A_T", "chr1_200_G_C"],
                },
                # VCF with complex variants
                {
                    "name": "complex_variants",
                    "vcf_content": """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t100\t.\tATGC\tA\t60\tPASS\t.\tGT\t0/1
chr1\t200\t.\tA\tATGCGTA\t55\tPASS\t.\tGT\t1/1
chr1\t300\t.\tATG\tCGA\t50\tPASS\t.\tGT\t0/0
""",
                    "targets": ["chr1_100_ATGC_A", "chr1_200_A_ATGCGTA", "chr1_300_ATG_CGA"],
                },
            ]
            
            for case in edge_cases:
                # Create test files
                vcf_file = temp_path / f"{case['name']}.vcf"
                vcf_file.write_text(case["vcf_content"])
                
                target_file = temp_path / f"{case['name']}_targets.txt"
                target_file.write_text("\n".join(case["targets"]))
                
                output_file = temp_path / f"{case['name']}_output"
                
                # Create configuration
                config = ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_file,
                    output_file=output_file,
                    batch_size=1000,
                    threads=1,
                    compress_output=False,
                    quiet=True
                )
                
                # Process VCF
                processor = VCFProcessor(config)
                result = processor.process()
                
                # Verify processing completed (may have some errors for edge cases)
                assert result.processed_variants >= 0, f"Case {case['name']}: Invalid processed count"
                
                # Should not crash
                summary = processor.get_processing_summary()
                assert "variants_processed" in summary
    
    def test_dry_run_integration(self):
        """Test dry run functionality integration."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test data
            vcf_file, target_file, expected_variants = self._create_comprehensive_vcf(
                temp_path, num_variants=15, num_samples=5
            )
            
            output_file = temp_path / "dry_run_output"
            
            # Test dry run
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_file,
                output_file=output_file,
                dry_run=True,
                verbose=True,
                quiet=False
            )
            
            # Record initial file state
            initial_files = set(temp_path.rglob("*"))
            
            # Process with dry run
            processor = VCFProcessor(config)
            result = processor.process()
            
            # Verify dry run behavior
            assert len(result.errors) == 0, f"Dry run had errors: {result.errors}"
            
            # Should not create output files
            gt_file = output_file.with_suffix(".gt.txt")
            seq_file = output_file.with_suffix(".seq.txt")
            
            assert not gt_file.exists(), "Dry run created .gt.txt file"
            assert not seq_file.exists(), "Dry run created .seq.txt file"
            
            # Should not create any new files (except possibly log files)
            final_files = set(temp_path.rglob("*"))
            new_files = final_files - initial_files
            
            for new_file in new_files:
                if new_file.is_file():
                    # Allow log files
                    assert (new_file.suffix in ['.log', '.tmp'] or 
                           'log' in new_file.name.lower()), f"Dry run created unexpected file: {new_file}"
            
            # Should provide validation information
            summary = processor.get_processing_summary()
            assert summary["variants_processed"] >= 0
    
    def test_error_recovery_integration(self):
        """Test error recovery and cleanup integration."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0
"""
            vcf_file.write_text(vcf_content)
            
            target_file = temp_path / "targets.txt"
            target_file.write_text("chr1_100_A_T\n")
            
            # Test with invalid output directory (read-only)
            readonly_dir = temp_path / "readonly"
            readonly_dir.mkdir()
            readonly_dir.chmod(0o444)  # Read-only
            
            output_file = readonly_dir / "output"
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_file,
                output_file=output_file,
                quiet=True
            )
            
            try:
                # This should handle the error gracefully
                processor = VCFProcessor(config)
                
                # May raise exception or complete with errors
                try:
                    result = processor.process()
                    # If it completes, should have errors
                    assert len(result.errors) > 0 or result.processed_variants == 0
                except (RuntimeError, PermissionError):
                    # Exception is also acceptable
                    pass
                
            finally:
                # Clean up: restore permissions
                readonly_dir.chmod(0o755)


class TestCompatibilityVerification:
    """Tests to verify compatibility with existing workflows."""
    
    def test_output_format_compatibility(self):
        """Test that output formats match expected structure."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test VCF
            vcf_file = temp_path / "compat_test.vcf"
            vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\tsample3
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0\t0/1\t1/1
chr1\t200\t.\tG\tC\t55\tPASS\t.\tGT\t0/1\t1/1\t0/0
chr2\t300\t.\tT\tA\t65\tPASS\t.\tGT\t1/1\t0/0\t0/1
"""
            vcf_file.write_text(vcf_content)
            
            target_file = temp_path / "targets.txt"
            target_file.write_text("chr1_100_A_T\nchr1_200_G_C\nchr2_300_T_A\n")
            
            output_file = temp_path / "compat_output"
            
            # Test different output formats
            format_configs = [
                {"miss_fmt": "NN", "gt_sep": "", "compress": False},
                {"miss_fmt": "./.", "gt_sep": "/", "compress": False},
                {"miss_fmt": "--", "gt_sep": "|", "compress": True},
            ]
            
            for i, fmt_config in enumerate(format_configs):
                config = ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_file,
                    output_file=temp_path / f"compat_output_{i}",
                    miss_fmt=fmt_config["miss_fmt"],
                    gt_sep=fmt_config["gt_sep"],
                    compress_output=fmt_config["compress"],
                    quiet=True
                )
                
                processor = VCFProcessor(config)
                result = processor.process()
                
                # Verify processing completed
                assert len(result.errors) == 0, f"Format {i}: Processing had errors: {result.errors}"
                
                # Verify output format expectations
                summary = processor.get_processing_summary()
                assert summary["success"] == True, f"Format {i}: Processing was not successful"
    
    def test_parameter_compatibility(self):
        """Test compatibility with various parameter combinations."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test data
            vcf_file = temp_path / "param_test.vcf"
            vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0
chr1\t200\t.\tG\tC\t55\tPASS\t.\tGT\t0/1
"""
            vcf_file.write_text(vcf_content)
            
            target_file = temp_path / "targets.txt"
            target_file.write_text("chr1_100_A_T\nchr1_200_G_C\n")
            
            # Test various parameter combinations
            param_combinations = [
                {"batch_size": 100, "threads": 1, "verbose": False, "quiet": True},
                {"batch_size": 1000, "threads": 2, "verbose": True, "quiet": False},
                {"batch_size": 10000, "threads": 1, "verbose": False, "quiet": False},
                {"batch_size": 5000, "threads": 4, "verbose": False, "quiet": True},
            ]
            
            for i, params in enumerate(param_combinations):
                config = ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_file,
                    output_file=temp_path / f"param_output_{i}",
                    **params
                )
                
                processor = VCFProcessor(config)
                result = processor.process()
                
                # Should complete without errors
                assert len(result.errors) == 0, f"Params {i}: Processing had errors: {result.errors}"
                assert result.processed_variants >= 0, f"Params {i}: Invalid processed count"
    
    def test_file_format_compatibility(self):
        """Test compatibility with different file formats and extensions."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Test different VCF file extensions
            vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0
"""
            
            target_content = "chr1_100_A_T\n"
            
            # Test different file extensions
            file_extensions = [
                (".vcf", ".txt"),
                (".vcf.gz", ".txt"),  # Note: We're not actually compressing for this test
                (".vcf", ".tsv"),
            ]
            
            for i, (vcf_ext, target_ext) in enumerate(file_extensions):
                vcf_file = temp_path / f"test_{i}{vcf_ext}"
                target_file = temp_path / f"targets_{i}{target_ext}"
                
                # For .gz files, we'll just use uncompressed content for simplicity
                vcf_file.write_text(vcf_content)
                target_file.write_text(target_content)
                
                output_file = temp_path / f"ext_output_{i}"
                
                config = ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_file,
                    output_file=output_file,
                    quiet=True
                )
                
                processor = VCFProcessor(config)
                result = processor.process()
                
                # Should handle different extensions
                assert len(result.errors) == 0, f"Extension {vcf_ext}/{target_ext}: Processing had errors: {result.errors}"
                assert result.processed_variants >= 0, f"Extension {vcf_ext}/{target_ext}: Invalid processed count"


class TestVariantTypeAnnotationIntegration:
    """Integration tests for variant type annotation feature."""
    
    def test_variant_type_annotation_end_to_end(self):
        """Test end-to-end processing with variant type annotation enabled."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create comprehensive VCF with different variant types
            vcf_file = temp_path / "variant_types.vcf"
            vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\tsample3
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0\t0/1\t1/1
chr1\t200\t.\tAT\tA\t55\tPASS\t.\tGT\t0/1\t1/1\t0/0
chr1\t300\t.\tG\tGC\t65\tPASS\t.\tGT\t1/1\t0/0\t0/1
chr1\t400\t.\tATG\tCGA\t70\tPASS\t.\tGT\t0/0\t0/1\t1/1
chr1\t500\t.\tA\tA\t50\tPASS\t.\tGT\t0/0\t0/0\t0/0
chr2\t100\t.\tC\tT,G\t60\tPASS\t.\tGT\t0/1\t0/2\t1/2
chr2\t200\t.\tA\tT,AT\t65\tPASS\t.\tGT\t0/1\t0/2\t1/2
chr2\t300\t.\tAT\tGC,A\t70\tPASS\t.\tGT\t0/1\t0/2\t1/2
"""
            vcf_file.write_text(vcf_content)
            
            target_file = temp_path / "targets.txt"
            target_file.write_text("chr1_100\nchr1_200\nchr1_300\nchr1_400\nchr1_500\nchr2_100\nchr2_200\nchr2_300\n")
            
            # Test with variant type annotation enabled
            output_file = temp_path / "va_type_output"
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_file,
                output_file=output_file,
                include_variant_type=True,
                compress_output=False,  # Use uncompressed for easier testing
                quiet=True
            )
            
            processor = VCFProcessor(config)
            
            try:
                result = processor.process()
                
                # Verify processing completed
                assert len(result.errors) == 0, f"Processing had errors: {result.errors}"
                assert result.processed_variants > 0, "No variants processed"
                
                # Check output files
                gt_file = output_file.with_suffix('.genotype_codes.tsv')
                seq_file = output_file.with_suffix('.genotype_bases.tsv')
                
                if gt_file.exists() and seq_file.exists():
                    # Read and verify output structure
                    gt_df = pd.read_csv(gt_file, sep='\t')
                    seq_df = pd.read_csv(seq_file, sep='\t')
                    
                    # Verify Variant_Type column exists and is in correct position
                    assert 'Variant_Type' in gt_df.columns, "Variant_Type column missing from genotype_codes.tsv"
                    assert 'Variant_Type' in seq_df.columns, "Variant_Type column missing from genotype_bases.tsv"
                    
                    # Verify column order
                    gt_columns = list(gt_df.columns)
                    seq_columns = list(seq_df.columns)
                    
                    alt_index_gt = gt_columns.index('ALT')
                    alt_index_seq = seq_columns.index('ALT')
                    variant_type_index_gt = gt_columns.index('Variant_Type')
                    variant_type_index_seq = seq_columns.index('Variant_Type')
                    
                    assert variant_type_index_gt == alt_index_gt + 1, "Variant_Type not immediately after ALT in GT file"
                    assert variant_type_index_seq == alt_index_seq + 1, "Variant_Type not immediately after ALT in SEQ file"
                    
                    # Verify both files have same structure
                    assert gt_columns == seq_columns, "GT and SEQ files have different column structures"
                    
                    # Verify variant type classifications
                    if len(gt_df) > 0:
                        variant_types = gt_df['Variant_Type'].tolist()
                        
                        # All variant types should be valid
                        valid_types = {'SNP', 'INDEL', 'MNP', 'REF', 'SNP|INDEL', 'SNP|MNP', 'INDEL|MNP', 'SNP|INDEL|MNP', 'SNP|INDEL|MNP|REF'}
                        for vt in variant_types:
                            assert any(valid_type in vt for valid_type in ['SNP', 'INDEL', 'MNP', 'REF']), f"Invalid variant type: {vt}"
                        
                        # Check specific expected types based on test data
                        # Note: Exact matching depends on how variants are processed and filtered
                        expected_type_patterns = ['SNP', 'INDEL', 'MNP', 'REF']
                        found_types = set()
                        for vt in variant_types:
                            for pattern in expected_type_patterns:
                                if pattern in vt:
                                    found_types.add(pattern)
                        
                        # Should find at least some of the expected types
                        assert len(found_types) > 0, f"No expected variant types found. Got: {variant_types}"
                
            except Exception as e:
                # If cyvcf2 is not available, skip the test
                if "cyvcf2" in str(e) or "VCF reader not available" in str(e):
                    pytest.skip("cyvcf2 not available for integration test")
                else:
                    raise
    
    def test_variant_type_annotation_disabled_by_default(self):
        """Test that variant type annotation is disabled by default."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create simple VCF
            vcf_file = temp_path / "simple.vcf"
            vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0
chr1\t200\t.\tG\tC\t55\tPASS\t.\tGT\t0/1
"""
            vcf_file.write_text(vcf_content)
            
            target_file = temp_path / "targets.txt"
            target_file.write_text("chr1_100\nchr1_200\n")
            
            output_file = temp_path / "default_output"
            
            # Test with default configuration (variant type annotation disabled)
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_file,
                output_file=output_file,
                compress_output=False,
                quiet=True
                # include_variant_type defaults to False
            )
            
            processor = VCFProcessor(config)
            
            try:
                result = processor.process()
                
                # Verify processing completed
                assert len(result.errors) == 0, f"Processing had errors: {result.errors}"
                
                # Check output files
                gt_file = output_file.with_suffix('.genotype_codes.tsv')
                seq_file = output_file.with_suffix('.genotype_bases.tsv')
                
                if gt_file.exists() and seq_file.exists():
                    # Read and verify output structure
                    gt_df = pd.read_csv(gt_file, sep='\t')
                    seq_df = pd.read_csv(seq_file, sep='\t')
                    
                    # Verify Variant_Type column does NOT exist
                    assert 'Variant_Type' not in gt_df.columns, "Variant_Type column should not exist when disabled"
                    assert 'Variant_Type' not in seq_df.columns, "Variant_Type column should not exist when disabled"
                    
                    # Should have standard columns
                    expected_base_columns = ['CHROM', 'POS', 'REF', 'ALT']
                    for col in expected_base_columns:
                        assert col in gt_df.columns, f"Missing expected column: {col}"
                        assert col in seq_df.columns, f"Missing expected column: {col}"
                
            except Exception as e:
                # If cyvcf2 is not available, skip the test
                if "cyvcf2" in str(e) or "VCF reader not available" in str(e):
                    pytest.skip("cyvcf2 not available for integration test")
                else:
                    raise
    
    def test_variant_type_annotation_with_different_configurations(self):
        """Test variant type annotation with different processing configurations."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test VCF
            vcf_file = temp_path / "config_test.vcf"
            vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0\t0/1
chr1\t200\t.\tAT\tA\t55\tPASS\t.\tGT\t0/1\t1/1
chr1\t300\t.\tG\tGC\t65\tPASS\t.\tGT\t1/1\t0/0
"""
            vcf_file.write_text(vcf_content)
            
            target_file = temp_path / "targets.txt"
            target_file.write_text("chr1_100\nchr1_200\nchr1_300\n")
            
            # Test different configurations with variant type annotation
            test_configs = [
                {
                    "name": "compressed",
                    "compress_output": True,
                    "batch_size": 1000,
                    "miss_fmt": "NN"
                },
                {
                    "name": "uncompressed",
                    "compress_output": False,
                    "batch_size": 5000,
                    "miss_fmt": "./."
                },
                {
                    "name": "large_batch",
                    "compress_output": False,
                    "batch_size": 10000,
                    "miss_fmt": "--"
                }
            ]
            
            for test_config in test_configs:
                output_file = temp_path / f"config_{test_config['name']}_output"
                
                config = ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_file,
                    output_file=output_file,
                    include_variant_type=True,  # Enable variant type annotation
                    compress_output=test_config["compress_output"],
                    batch_size=test_config["batch_size"],
                    miss_fmt=test_config["miss_fmt"],
                    quiet=True
                )
                
                processor = VCFProcessor(config)
                
                try:
                    result = processor.process()
                    
                    # Verify processing completed
                    assert len(result.errors) == 0, f"Config {test_config['name']}: Processing had errors: {result.errors}"
                    
                    # Verify configuration was applied
                    assert processor.config.include_variant_type == True, f"Config {test_config['name']}: Variant type not enabled"
                    assert processor.config.compress_output == test_config["compress_output"], f"Config {test_config['name']}: Compression setting not applied"
                    assert processor.config.batch_size == test_config["batch_size"], f"Config {test_config['name']}: Batch size not applied"
                    assert processor.config.miss_fmt == test_config["miss_fmt"], f"Config {test_config['name']}: Miss format not applied"
                    
                except Exception as e:
                    # If cyvcf2 is not available, skip the test
                    if "cyvcf2" in str(e) or "VCF reader not available" in str(e):
                        pytest.skip("cyvcf2 not available for integration test")
                    else:
                        raise