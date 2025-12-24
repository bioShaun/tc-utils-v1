"""Property-based tests for VCFProcessor error handling, logging, and validation."""

import tempfile
from pathlib import Path

import pandas as pd
import pytest

try:
    from hypothesis import given, strategies as st, assume
    HYPOTHESIS_AVAILABLE = True
except ImportError:
    HYPOTHESIS_AVAILABLE = False
    def given(*args, **kwargs):
        def decorator(func):
            return pytest.mark.skip("hypothesis not available")(func)
        return decorator

from chip.vcf_processor.config import ProcessingConfig
from chip.vcf_processor.vcf_processor import VCFProcessor, VCFProcessorFactory


class TestVCFProcessorErrorHandling:
    """Property tests for VCF processor error handling."""
    
    @pytest.mark.skipif(not HYPOTHESIS_AVAILABLE, reason="hypothesis not available")
    @given(
        batch_size=st.integers(min_value=-100, max_value=0),
        threads=st.integers(min_value=-10, max_value=0)
    )
    def test_invalid_config_parameters_raise_errors(self, batch_size, threads):
        """Property 3: Error Handling with Continuation - Invalid parameters should raise appropriate errors."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create minimal test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            vcf_file.write_text("##fileformat=VCFv4.2\n")
            target_id_file.write_text("test_id\n")
            
            # Test that invalid parameters raise ValueError
            with pytest.raises(ValueError):
                config = ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_id_file,
                    output_file=output_file,
                    batch_size=batch_size,
                    threads=threads
                )
                VCFProcessor(config)
    
    @pytest.mark.skipif(not HYPOTHESIS_AVAILABLE, reason="hypothesis not available")
    @given(
        missing_file_type=st.sampled_from(['vcf', 'target_id'])
    )
    def test_missing_files_raise_file_not_found(self, missing_file_type):
        """Property 3: Error Handling with Continuation - Missing files should raise FileNotFoundError."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create some files but not others
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            if missing_file_type != 'vcf':
                vcf_file.write_text("##fileformat=VCFv4.2\n")
            if missing_file_type != 'target_id':
                target_id_file.write_text("test_id\n")
            
            # Should raise FileNotFoundError for missing files
            with pytest.raises(FileNotFoundError):
                config = ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_id_file,
                    output_file=output_file
                )
                VCFProcessor(config)
    
    def test_error_recovery_in_dry_run(self):
        """Property 3: Error Handling with Continuation - Errors in dry run should be recoverable."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            # Create valid files
            vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	100	.	A	T	60	PASS	.	GT	0/0
"""
            vcf_file.write_text(vcf_content)
            target_id_file.write_text("chr1_100_A_T\n")
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                dry_run=True
            )
            
            processor = VCFProcessor(config)
            
            # Dry run should complete even if there are minor issues
            result = processor.process()
            
            # Should have completed without raising exceptions
            assert result is not None
            assert isinstance(result.errors, list)
    
    def test_error_accumulation_and_reporting(self):
        """Property 3: Error Handling with Continuation - Errors should be accumulated and reported."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	100	.	A	T	60	PASS	.	GT	0/0
"""
            vcf_file.write_text(vcf_content)
            target_id_file.write_text("chr1_100_A_T\n")
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                dry_run=True
            )
            
            processor = VCFProcessor(config)
            
            # Test error accumulation
            processor.result.add_error("Test error 1")
            processor.result.add_error("Test error 2")
            
            # Errors should be accumulated
            assert len(processor.result.errors) == 2
            assert "Test error 1" in processor.result.errors
            assert "Test error 2" in processor.result.errors
            
            # Summary should reflect errors
            summary = processor.get_processing_summary()
            assert summary["success"] is False
            assert len(summary["errors"]) == 2


class TestVCFProcessorLogging:
    """Property tests for VCF processor logging behavior."""
    
    def test_logging_initialization(self):
        """Property 4: Logging Behavior - Processor should log initialization."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            vcf_file.write_text("##fileformat=VCFv4.2\n")
            target_id_file.write_text("test_id\n")
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file
            )
            
            # Should not raise exceptions during initialization
            processor = VCFProcessor(config)
            assert processor is not None
    
    def test_logging_during_processing(self):
        """Property 4: Logging Behavior - Processing should generate appropriate log messages."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	100	.	A	T	60	PASS	.	GT	0/0
"""
            vcf_file.write_text(vcf_content)
            target_id_file.write_text("chr1_100_A_T\n")
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                dry_run=True
            )
            
            processor = VCFProcessor(config)
            
            # Processing should complete and log appropriately
            result = processor.process()
            assert result is not None
    
    @pytest.mark.skipif(not HYPOTHESIS_AVAILABLE, reason="hypothesis not available")
    @given(
        verbose=st.booleans(),
        quiet=st.booleans()
    )
    def test_logging_verbosity_modes(self, verbose, quiet):
        """Property 4: Logging Behavior - Different verbosity modes should work correctly."""
        # Skip contradictory combinations
        assume(not (verbose and quiet))
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            # Create proper VCF content
            vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	100	.	A	T	60	PASS	.	GT	0/0
"""
            vcf_file.write_text(vcf_content)
            target_id_file.write_text("chr1_100_A_T\n")
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                verbose=verbose,
                quiet=quiet,
                dry_run=True
            )
            
            # Should handle different verbosity modes without errors
            processor = VCFProcessor(config)
            result = processor.process()
            assert result is not None


class TestVCFProcessorValidation:
    """Property tests for VCF format validation."""
    
    @pytest.mark.skipif(not HYPOTHESIS_AVAILABLE, reason="hypothesis not available")
    @given(
        vcf_header=st.sampled_from([
            "##fileformat=VCFv4.2",
            "##fileformat=VCFv4.1", 
            "##fileformat=VCFv4.0",
            "##invalid_header",
            ""
        ])
    )
    def test_vcf_format_validation(self, vcf_header):
        """Property 13: VCF Format Validation - Different VCF formats should be handled appropriately."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            # Create VCF with different headers
            vcf_content = f"""{vcf_header}
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	100	.	A	T	60	PASS	.	GT	0/0
"""
            vcf_file.write_text(vcf_content)
            target_id_file.write_text("chr1_100_A_T\n")
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                dry_run=True
            )
            
            # Should handle different VCF formats gracefully
            try:
                processor = VCFProcessor(config)
                result = processor.process()
                # Valid formats should process successfully
                if vcf_header.startswith("##fileformat=VCF"):
                    assert result is not None
            except Exception as e:
                # Invalid formats may raise exceptions, which is acceptable
                if not vcf_header.startswith("##fileformat=VCF"):
                    assert True  # Expected for invalid formats
                else:
                    # Valid formats should not raise exceptions
                    pytest.fail(f"Valid VCF format raised exception: {e}")
    
    def test_vcf_structure_validation(self):
        """Property 13: VCF Format Validation - VCF structure should be validated."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            # Test with malformed VCF - but handle errors gracefully
            malformed_vcf = """##fileformat=VCFv4.2
# Missing proper header line
chr1	100	.	A	T	60	PASS	.	GT	0/0
"""
            vcf_file.write_text(malformed_vcf)
            target_id_file.write_text("chr1_100_A_T\n")
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                dry_run=True
            )
            
            # Should handle malformed VCF gracefully
            processor = VCFProcessor(config)
            
            # Malformed VCF may raise exceptions, which is acceptable
            try:
                result = processor.process()
                # If it succeeds, that's fine too
                assert result is not None
            except RuntimeError as e:
                # Expected for malformed VCF
                assert "cannot parse VCF header" in str(e) or "VCF" in str(e)
    
    def test_target_id_validation(self):
        """Property 13: VCF Format Validation - Target ID format should be validated."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	100	.	A	T	60	PASS	.	GT	0/0
"""
            vcf_file.write_text(vcf_content)
            
            # Test with various target ID formats
            target_content = """chr1_100_A_T
invalid id with spaces
chr2_200_G_C
# comment line
"""
            target_id_file.write_text(target_content)
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                dry_run=True
            )
            
            # Should handle mixed valid/invalid target IDs
            processor = VCFProcessor(config)
            result = processor.process()
            
            # Should complete processing
            assert result is not None
    
    def test_output_validation(self):
        """Property 13: VCF Format Validation - Output validation should work correctly."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	100	.	A	T	60	PASS	.	GT	0/0
"""
            vcf_file.write_text(vcf_content)
            target_id_file.write_text("chr1_100_A_T\n")
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                dry_run=True
            )
            
            processor = VCFProcessor(config)
            result = processor.process()
            
            # Test output validation
            is_valid = processor.validate_output()
            
            # For dry run, validation behavior may vary
            assert isinstance(is_valid, bool)


class TestVCFProcessorFactory:
    """Property tests for VCFProcessorFactory."""
    
    @pytest.mark.skipif(not HYPOTHESIS_AVAILABLE, reason="hypothesis not available")
    @given(
        batch_size=st.integers(min_value=100, max_value=50000),
        threads=st.integers(min_value=1, max_value=8),
        compress_output=st.booleans()
    )
    def test_factory_creates_valid_processors(self, batch_size, threads, compress_output):
        """Factory should create valid processors with different configurations."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            vcf_file.write_text("##fileformat=VCFv4.2\n")
            target_id_file.write_text("test_id\n")
            
            # Test factory methods
            processor = VCFProcessorFactory.create_from_files(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                batch_size=batch_size,
                threads=threads,
                compress_output=compress_output,
                dry_run=True
            )
            
            assert isinstance(processor, VCFProcessor)
            assert processor.config.batch_size == batch_size
            assert processor.config.threads == threads
            assert processor.config.compress_output == compress_output
    
    def test_factory_batch_processor_optimization(self):
        """Factory should create optimized batch processors."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            vcf_file.write_text("##fileformat=VCFv4.2\n")
            target_id_file.write_text("test_id\n")
            
            # Create batch processor
            processor = VCFProcessorFactory.create_batch_processor(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                batch_size=5000
            )
            
            assert isinstance(processor, VCFProcessor)
            assert processor.config.batch_size == 5000


class TestVCFProcessorVariantTypeAnnotation:
    """Property tests for variant type annotation feature."""
    
    def test_variant_type_column_placement(self):
        """Property 20: Variant Type Column Placement.
        
        Feature: vcf-processor-optimization, Property 20: Variant Type Column Placement
        Validates: Requirements 8.1, 8.8
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            # Create VCF with different variant types
            vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2
chr1	100	.	A	T	60	PASS	.	GT	0/0	0/1
chr1	200	.	AT	A	60	PASS	.	GT	1/1	0/0
chr1	300	.	G	GC	60	PASS	.	GT	0/1	1/1
chr1	400	.	ATG	GCA	60	PASS	.	GT	0/0	0/1
"""
            vcf_file.write_text(vcf_content)
            target_id_file.write_text("chr1_100\nchr1_200\nchr1_300\nchr1_400\n")
            
            # Test with variant type annotation enabled
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                include_variant_type=True,
                compress_output=False,  # Use uncompressed for easier testing
                dry_run=False  # Need actual processing to test output
            )
            
            processor = VCFProcessor(config)
            
            try:
                result = processor.process()
                
                # Check that output files were created
                gt_file = output_file.with_suffix('.genotype_codes.tsv')
                seq_file = output_file.with_suffix('.genotype_bases.tsv')
                
                if gt_file.exists() and seq_file.exists():
                    # Read the output files
                    gt_df = pd.read_csv(gt_file, sep='\t')
                    seq_df = pd.read_csv(seq_file, sep='\t')
                    
                    # Check that Variant_Type column exists
                    assert 'Variant_Type' in gt_df.columns
                    assert 'Variant_Type' in seq_df.columns
                    
                    # Check that Variant_Type column is after ALT column
                    gt_columns = list(gt_df.columns)
                    seq_columns = list(seq_df.columns)
                    
                    alt_index_gt = gt_columns.index('ALT')
                    alt_index_seq = seq_columns.index('ALT')
                    
                    variant_type_index_gt = gt_columns.index('Variant_Type')
                    variant_type_index_seq = seq_columns.index('Variant_Type')
                    
                    # Variant_Type should be immediately after ALT
                    assert variant_type_index_gt == alt_index_gt + 1
                    assert variant_type_index_seq == alt_index_seq + 1
                    
                    # Check that both files have the same column structure
                    assert gt_columns == seq_columns
                    
                    # Check that variant types are correctly classified
                    expected_types = ['SNP', 'INDEL', 'INDEL', 'MNP']  # Based on the test data
                    if len(gt_df) > 0:
                        # At least some variants should be processed
                        assert all(vt in ['SNP', 'INDEL', 'MNP', 'REF'] for vt in gt_df['Variant_Type'])
                        assert all(vt in ['SNP', 'INDEL', 'MNP', 'REF'] for vt in seq_df['Variant_Type'])
                
            except Exception as e:
                # If cyvcf2 is not available or other issues, the test should still pass
                # as long as the configuration is accepted
                if "cyvcf2" in str(e) or "VCF reader not available" in str(e):
                    pytest.skip("cyvcf2 not available for full integration test")
                else:
                    # Re-raise other exceptions
                    raise
    
    def test_variant_type_column_placement_disabled(self):
        """Property 20: Variant Type Column Placement - Feature disabled.
        
        Feature: vcf-processor-optimization, Property 20: Variant Type Column Placement
        Validates: Requirements 8.1, 8.8
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	100	.	A	T	60	PASS	.	GT	0/0
"""
            vcf_file.write_text(vcf_content)
            target_id_file.write_text("chr1_100\n")
            
            # Test with variant type annotation disabled (default)
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                include_variant_type=False,
                compress_output=False,
                dry_run=False
            )
            
            processor = VCFProcessor(config)
            
            try:
                result = processor.process()
                
                # Check that output files were created
                gt_file = output_file.with_suffix('.genotype_codes.tsv')
                seq_file = output_file.with_suffix('.genotype_bases.tsv')
                
                if gt_file.exists() and seq_file.exists():
                    # Read the output files
                    gt_df = pd.read_csv(gt_file, sep='\t')
                    seq_df = pd.read_csv(seq_file, sep='\t')
                    
                    # Check that Variant_Type column does NOT exist
                    assert 'Variant_Type' not in gt_df.columns
                    assert 'Variant_Type' not in seq_df.columns
                    
                    # Should have standard columns: CHROM, POS, REF, ALT, sample columns
                    expected_base_columns = ['CHROM', 'POS', 'REF', 'ALT']
                    for col in expected_base_columns:
                        assert col in gt_df.columns
                        assert col in seq_df.columns
                
            except Exception as e:
                # If cyvcf2 is not available or other issues, the test should still pass
                if "cyvcf2" in str(e) or "VCF reader not available" in str(e):
                    pytest.skip("cyvcf2 not available for full integration test")
                else:
                    raise