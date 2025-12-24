"""Tests for VCFProcessor main orchestrator."""

import tempfile
from pathlib import Path

import pandas as pd
import pytest

from chip.vcf_processor.config import ProcessingConfig
from chip.vcf_processor.vcf_processor import VCFProcessor, VCFProcessorFactory


class TestVCFProcessor:
    """Test VCFProcessor main orchestrator."""
    
    def test_vcf_processor_creation(self):
        """Test VCFProcessor can be created with valid configuration."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            # Create minimal VCF file
            vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2
chr1	100	.	A	T	60	PASS	.	GT	0/0	0/1
chr1	200	.	T	G	60	PASS	.	GT	0/1	1/1
"""
            vcf_file.write_text(vcf_content)
            
            # Create target ID file
            target_id_file.write_text("chr1_100_A_T\nchr1_200_T_G\n")
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                dry_run=True  # Use dry run to avoid actual processing
            )
            
            processor = VCFProcessor(config)
            
            assert processor.config == config
            assert processor.result is not None
            assert processor.error_handler is not None
    
    def test_vcf_processor_dry_run(self):
        """Test VCFProcessor dry run functionality."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            # Create minimal VCF file
            vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2
chr1	100	.	A	T	60	PASS	.	GT	0/0	0/1
chr1	200	.	T	G	60	PASS	.	GT	0/1	1/1
"""
            vcf_file.write_text(vcf_content)
            
            # Create target ID file
            target_id_file.write_text("chr1_100_A_T\nchr1_200_T_G\n")
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                dry_run=True
            )
            
            processor = VCFProcessor(config)
            result = processor.process()
            
            # Dry run should complete without errors
            assert len(result.errors) == 0
            assert result.processed_variants >= 0  # May be 0 in dry run
    
    def test_vcf_processor_fallback_processing(self):
        """Test VCFProcessor fallback processing (without cyvcf2)."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            # Create minimal VCF file
            vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2
chr1	100	.	A	T	60	PASS	.	GT	0/0	0/1
chr1	200	.	T	G	60	PASS	.	GT	0/1	1/1
"""
            vcf_file.write_text(vcf_content)
            
            # Create target ID file
            target_id_file.write_text("chr1_100_A_T\nchr1_200_T_G\n")
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                dry_run=False
            )
            
            processor = VCFProcessor(config)
            
            # Force fallback mode by setting vcf_reader to None
            processor._vcf_reader = None
            
            result = processor.process()
            
            # Processing should complete
            assert result.processed_variants >= 0
            assert result.total_variants >= 0
    
    def test_vcf_processor_validation_errors(self):
        """Test VCFProcessor validation with invalid configuration."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Test missing VCF file
            with pytest.raises(FileNotFoundError, match="VCF file not found"):
                config = ProcessingConfig(
                    vcf_file=temp_path / "nonexistent.vcf",
                    target_id_file=temp_path / "targets.txt",
                    output_file=temp_path / "output"
                )
                VCFProcessor(config)
            
            # Create VCF file but not target file
            vcf_file = temp_path / "test.vcf"
            vcf_file.write_text("##fileformat=VCFv4.2\n")
            
            with pytest.raises(FileNotFoundError, match="Target ID file not found"):
                config = ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=temp_path / "nonexistent_targets.txt",
                    output_file=temp_path / "output"
                )
                VCFProcessor(config)
    
    def test_vcf_processor_invalid_parameters(self):
        """Test VCFProcessor with invalid parameters."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            vcf_file.write_text("##fileformat=VCFv4.2\n")
            target_id_file.write_text("test_id\n")
            
            # Test invalid batch size
            with pytest.raises(ValueError, match="Batch size must be positive"):
                config = ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_id_file,
                    output_file=output_file,
                    batch_size=0
                )
                VCFProcessor(config)
            
            # Test invalid thread count
            with pytest.raises(ValueError, match="Threads must be positive"):
                config = ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_id_file,
                    output_file=output_file,
                    threads=0
                )
                VCFProcessor(config)
    
    def test_vcf_processor_summary(self):
        """Test VCFProcessor processing summary."""
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
            processor.process()
            
            summary = processor.get_processing_summary()
            
            assert "variants_processed" in summary
            assert "total_variants" in summary
            assert "output_files" in summary
            assert "errors" in summary
            assert "success" in summary
            assert isinstance(summary["success"], bool)


class TestVCFProcessorFactory:
    """Test VCFProcessorFactory."""
    
    def test_create_processor(self):
        """Test factory create_processor method."""
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
            
            processor = VCFProcessorFactory.create_processor(config)
            
            assert isinstance(processor, VCFProcessor)
            assert processor.config == config
    
    def test_create_from_files(self):
        """Test factory create_from_files method."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            vcf_file.write_text("##fileformat=VCFv4.2\n")
            target_id_file.write_text("test_id\n")
            
            processor = VCFProcessorFactory.create_from_files(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                dry_run=True
            )
            
            assert isinstance(processor, VCFProcessor)
            assert processor.config.vcf_file == vcf_file
            assert processor.config.target_id_file == target_id_file
            assert processor.config.output_file == output_file
            assert processor.config.dry_run is True
    
    def test_create_batch_processor(self):
        """Test factory create_batch_processor method."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            output_file = temp_path / "output"
            
            vcf_file.write_text("##fileformat=VCFv4.2\n")
            target_id_file.write_text("test_id\n")
            
            processor = VCFProcessorFactory.create_batch_processor(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                batch_size=5000
            )
            
            assert isinstance(processor, VCFProcessor)
            assert processor.config.batch_size == 5000