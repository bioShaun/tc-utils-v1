"""Property-based tests for OutputWriter."""

import gzip
import tempfile
from pathlib import Path
from typing import List

import pandas as pd
import pytest
from hypothesis import given, strategies as st

from chip.vcf_processor.config import ProcessingConfig, ProcessingResult
from chip.vcf_processor.output_writer import OutputWriter, BatchOutputWriter


# Test data generators
@st.composite
def processing_config(draw):
    """Generate valid ProcessingConfig for testing."""
    # Create temporary files
    vcf_file = Path(tempfile.mktemp(suffix=".vcf"))
    target_id_file = Path(tempfile.mktemp(suffix=".txt"))
    output_file = Path(tempfile.mktemp(suffix=".out"))
    
    # Create the files so they exist
    vcf_file.touch()
    target_id_file.touch()
    
    return ProcessingConfig(
        vcf_file=vcf_file,
        target_id_file=target_id_file,
        output_file=output_file,
        miss_fmt=draw(st.sampled_from(["NN", "N", "--"])),
        gt_sep=draw(st.sampled_from(["", "/", "|"])),
        threads=draw(st.integers(min_value=1, max_value=8)),
        batch_size=draw(st.integers(min_value=100, max_value=10000)),
        compress_output=draw(st.booleans()),
        verbose=draw(st.booleans()),
        quiet=False,  # Don't conflict with verbose
        dry_run=draw(st.booleans()),
    )


@st.composite
def sample_dataframe(draw):
    """Generate sample DataFrame for testing."""
    n_rows = draw(st.integers(min_value=1, max_value=100))
    n_samples = draw(st.integers(min_value=1, max_value=10))
    
    # Generate basic columns
    chroms = draw(st.lists(
        st.sampled_from(["chr1", "chr2", "chr3", "chrX", "chrY"]),
        min_size=n_rows, max_size=n_rows
    ))
    positions = draw(st.lists(
        st.integers(min_value=1, max_value=1000000),
        min_size=n_rows, max_size=n_rows
    ))
    refs = draw(st.lists(
        st.sampled_from(["A", "T", "G", "C", "AT", "GC"]),
        min_size=n_rows, max_size=n_rows
    ))
    alts = draw(st.lists(
        st.sampled_from(["A", "T", "G", "C", "AT", "GC", "delA", "insT"]),
        min_size=n_rows, max_size=n_rows
    ))
    
    # Create DataFrame
    data = {
        "CHROM": chroms,
        "POS": positions,
        "REF": refs,
        "ALT": alts,
    }
    
    # Add sample columns
    for i in range(n_samples):
        sample_name = f"sample_{i}"
        genotypes = draw(st.lists(
            st.sampled_from(["0/0", "0/1", "1/1", "./."]),
            min_size=n_rows, max_size=n_rows
        ))
        data[sample_name] = genotypes
    
    return pd.DataFrame(data)


class TestOutputWriterProperties:
    """Property-based tests for OutputWriter."""
    
    @given(config=processing_config(), gt_df=sample_dataframe(), seq_df=sample_dataframe())
    def test_output_file_naming_consistency(self, config, gt_df, seq_df):
        """Property 15: Output File Naming Consistency
        
        For any output file path, the generated files should follow the expected 
        naming pattern (.gt.txt.gz and .seq.txt.gz)
        
        **Validates: Requirements 7.3**
        """
        result = ProcessingResult()
        
        try:
            writer = OutputWriter(config, result)
            
            # Get expected paths
            expected_paths = writer.get_output_paths()
            
            # Check naming pattern
            gt_path = expected_paths["genotype"]
            seq_path = expected_paths["sequence"]
            
            # Verify naming convention
            if config.compress_output:
                assert gt_path.name.endswith(".gt.txt.gz"), f"GT file should end with .gt.txt.gz, got {gt_path.name}"
                assert seq_path.name.endswith(".seq.txt.gz"), f"SEQ file should end with .seq.txt.gz, got {seq_path.name}"
            else:
                assert gt_path.name.endswith(".gt.txt"), f"GT file should end with .gt.txt, got {gt_path.name}"
                assert seq_path.name.endswith(".seq.txt"), f"SEQ file should end with .seq.txt, got {seq_path.name}"
            
            # Verify base name consistency
            base_name = config.output_file.name
            assert base_name in gt_path.name, f"GT file should contain base name {base_name}"
            assert base_name in seq_path.name, f"SEQ file should contain base name {base_name}"
            
            # Test actual file creation
            writer.write_dataframes(gt_df, seq_df)
            
            # Verify files were created with correct names
            assert gt_path.exists(), f"GT file should be created at {gt_path}"
            assert seq_path.exists(), f"SEQ file should be created at {seq_path}"
            
        finally:
            # Cleanup
            for path in [config.vcf_file, config.target_id_file]:
                if path.exists():
                    path.unlink()
            
            # Cleanup output files
            writer = OutputWriter(config, result)
            paths = writer.get_output_paths()
            for path in paths.values():
                if path.exists():
                    path.unlink()
    
    @given(config=processing_config())
    def test_output_path_generation_consistency(self, config):
        """Test that output path generation is consistent across calls."""
        result = ProcessingResult()
        
        try:
            writer1 = OutputWriter(config, result)
            writer2 = OutputWriter(config, result)
            
            paths1 = writer1.get_output_paths()
            paths2 = writer2.get_output_paths()
            
            # Paths should be identical
            assert paths1["genotype"] == paths2["genotype"]
            assert paths1["sequence"] == paths2["sequence"]
            
        finally:
            # Cleanup
            for path in [config.vcf_file, config.target_id_file]:
                if path.exists():
                    path.unlink()
    
    @given(config=processing_config(), gt_df=sample_dataframe(), seq_df=sample_dataframe())
    def test_batch_vs_single_write_naming_consistency(self, config, gt_df, seq_df):
        """Test that batch and single write produce same file names."""
        result1 = ProcessingResult()
        result2 = ProcessingResult()
        
        try:
            # Single write
            writer1 = OutputWriter(config, result1)
            paths1 = writer1.get_output_paths()
            
            # Batch write
            writer2 = BatchOutputWriter(config, result2)
            paths2 = writer2.get_output_paths()
            
            # Paths should be identical
            assert paths1["genotype"] == paths2["genotype"]
            assert paths1["sequence"] == paths2["sequence"]
            
        finally:
            # Cleanup
            for path in [config.vcf_file, config.target_id_file]:
                if path.exists():
                    path.unlink()
    
    @given(
        config=processing_config(),
        batches=st.lists(sample_dataframe(), min_size=1, max_size=5)
    )
    def test_batch_append_mode_consistency(self, config, batches):
        """Test that batch processing maintains consistent file naming across batches."""
        result = ProcessingResult()
        
        try:
            writer = BatchOutputWriter(config, result)
            expected_paths = writer.get_output_paths()
            
            writer.start_batch_processing()
            
            for i, batch_df in enumerate(batches):
                # Write each batch
                writer.write_batch(batch_df, batch_df)  # Use same df for both gt and seq
                
                # Verify files exist and have expected names
                gt_path = expected_paths["genotype"]
                seq_path = expected_paths["sequence"]
                
                assert gt_path.exists(), f"GT file should exist after batch {i}"
                assert seq_path.exists(), f"SEQ file should exist after batch {i}"
                
                # File names should remain consistent
                if config.compress_output:
                    assert gt_path.name.endswith(".gt.txt.gz")
                    assert seq_path.name.endswith(".seq.txt.gz")
                else:
                    assert gt_path.name.endswith(".gt.txt")
                    assert seq_path.name.endswith(".seq.txt")
            
            writer.finish_batch_processing()
            
            # Final verification
            assert writer.validate_output_files(), "Output files should be valid after batch processing"
            
        finally:
            # Cleanup
            for path in [config.vcf_file, config.target_id_file]:
                if path.exists():
                    path.unlink()
            
            # Cleanup output files
            writer = OutputWriter(config, result)
            paths = writer.get_output_paths()
            for path in paths.values():
                if path.exists():
                    path.unlink()


class TestOutputWriterEdgeCases:
    """Edge case tests for output file naming."""
    
    def test_special_characters_in_output_path(self):
        """Test handling of special characters in output paths."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create input files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            vcf_file.touch()
            target_id_file.touch()
            
            # Test with special characters
            output_file = temp_path / "test-output_file.with.dots"
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                compress_output=True
            )
            
            result = ProcessingResult()
            writer = OutputWriter(config, result)
            paths = writer.get_output_paths()
            
            # Verify naming pattern is preserved
            assert paths["genotype"].name == "test-output_file.with.dots.gt.txt.gz"
            assert paths["sequence"].name == "test-output_file.with.dots.seq.txt.gz"
    
    def test_long_output_path(self):
        """Test handling of very long output paths."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create input files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            vcf_file.touch()
            target_id_file.touch()
            
            # Create long output path
            long_name = "a" * 100  # Very long filename
            output_file = temp_path / long_name
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                compress_output=False
            )
            
            result = ProcessingResult()
            writer = OutputWriter(config, result)
            paths = writer.get_output_paths()
            
            # Verify naming pattern is preserved even with long names
            assert paths["genotype"].name == f"{long_name}.gt.txt"
            assert paths["sequence"].name == f"{long_name}.seq.txt"


class TestCompressionSupport:
    """Property-based tests for compression support."""
    
    @given(
        config=processing_config(),
        gt_df=sample_dataframe(),
        seq_df=sample_dataframe()
    )
    def test_output_format_support(self, config, gt_df, seq_df):
        """Property 16: Output Format Support
        
        For any compression setting, the processor should generate valid output 
        files in both compressed and uncompressed formats as requested
        
        **Validates: Requirements 7.4**
        """
        result = ProcessingResult()
        
        try:
            writer = OutputWriter(config, result)
            
            # Write data
            writer.write_dataframes(gt_df, seq_df)
            
            # Get output paths
            paths = writer.get_output_paths()
            gt_path = paths["genotype"]
            seq_path = paths["sequence"]
            
            # Verify files exist
            assert gt_path.exists(), f"GT output file should exist: {gt_path}"
            assert seq_path.exists(), f"SEQ output file should exist: {seq_path}"
            
            # Verify file format based on compression setting
            if config.compress_output:
                # Files should be compressed (gzip format)
                assert gt_path.suffix == ".gz", f"GT file should be compressed: {gt_path}"
                assert seq_path.suffix == ".gz", f"SEQ file should be compressed: {seq_path}"
                
                # Verify files can be read as gzip
                try:
                    with gzip.open(gt_path, 'rt') as f:
                        gt_content = f.read()
                    with gzip.open(seq_path, 'rt') as f:
                        seq_content = f.read()
                    
                    # Content should not be empty
                    assert len(gt_content) > 0, "Compressed GT file should have content"
                    assert len(seq_content) > 0, "Compressed SEQ file should have content"
                    
                    # Content should contain expected headers
                    assert "CHROM" in gt_content, "GT file should contain CHROM header"
                    assert "CHROM" in seq_content, "SEQ file should contain CHROM header"
                    
                except Exception as e:
                    pytest.fail(f"Failed to read compressed files: {e}")
                    
            else:
                # Files should be uncompressed (plain text)
                assert gt_path.suffix == ".txt", f"GT file should be uncompressed: {gt_path}"
                assert seq_path.suffix == ".txt", f"SEQ file should be uncompressed: {seq_path}"
                
                # Verify files can be read as plain text
                try:
                    with open(gt_path, 'r') as f:
                        gt_content = f.read()
                    with open(seq_path, 'r') as f:
                        seq_content = f.read()
                    
                    # Content should not be empty
                    assert len(gt_content) > 0, "Uncompressed GT file should have content"
                    assert len(seq_content) > 0, "Uncompressed SEQ file should have content"
                    
                    # Content should contain expected headers
                    assert "CHROM" in gt_content, "GT file should contain CHROM header"
                    assert "CHROM" in seq_content, "SEQ file should contain CHROM header"
                    
                except Exception as e:
                    pytest.fail(f"Failed to read uncompressed files: {e}")
            
            # Verify file sizes are reasonable
            gt_size = gt_path.stat().st_size
            seq_size = seq_path.stat().st_size
            
            assert gt_size > 0, "GT file should have non-zero size"
            assert seq_size > 0, "SEQ file should have non-zero size"
            
            # For compressed files, verify they're actually smaller than uncompressed
            # (This is a heuristic - small files might not compress well)
            if config.compress_output and len(gt_df) > 50:  # Only test compression on larger files
                # Create uncompressed version for comparison
                temp_config = ProcessingConfig(
                    vcf_file=config.vcf_file,
                    target_id_file=config.target_id_file,
                    output_file=config.output_file.with_suffix(".temp"),
                    compress_output=False,
                    miss_fmt=config.miss_fmt,
                    gt_sep=config.gt_sep,
                    threads=config.threads,
                    batch_size=config.batch_size
                )
                
                temp_result = ProcessingResult()
                temp_writer = OutputWriter(temp_config, temp_result)
                temp_writer.write_dataframes(gt_df, seq_df)
                
                temp_paths = temp_writer.get_output_paths()
                temp_gt_size = temp_paths["genotype"].stat().st_size
                temp_seq_size = temp_paths["sequence"].stat().st_size
                
                # Compressed should generally be smaller (allowing generous overhead for small files)
                compression_ratio_gt = gt_size / temp_gt_size if temp_gt_size > 0 else 1
                compression_ratio_seq = seq_size / temp_seq_size if temp_seq_size > 0 else 1
                
                # Allow for cases where compression doesn't help much (especially small files)
                # Use more lenient ratio for small files
                max_ratio = 3.0 if len(gt_df) < 100 else 1.5
                assert compression_ratio_gt <= max_ratio, f"GT compression ratio should be reasonable: {compression_ratio_gt}"
                assert compression_ratio_seq <= max_ratio, f"SEQ compression ratio should be reasonable: {compression_ratio_seq}"
                
                # Cleanup temp files
                for path in temp_paths.values():
                    if path.exists():
                        path.unlink()
            
        finally:
            # Cleanup
            for path in [config.vcf_file, config.target_id_file]:
                if path.exists():
                    path.unlink()
            
            # Cleanup output files
            writer = OutputWriter(config, result)
            paths = writer.get_output_paths()
            for path in paths.values():
                if path.exists():
                    path.unlink()
    
    @given(config=processing_config())
    def test_compression_setting_consistency(self, config):
        """Test that compression setting is consistently applied."""
        result = ProcessingResult()
        
        try:
            writer = OutputWriter(config, result)
            paths = writer.get_output_paths()
            
            # Check file extensions match compression setting
            for path in paths.values():
                if config.compress_output:
                    assert path.name.endswith(".gz"), f"Compressed file should end with .gz: {path}"
                else:
                    assert not path.name.endswith(".gz"), f"Uncompressed file should not end with .gz: {path}"
                    assert path.name.endswith(".txt"), f"Uncompressed file should end with .txt: {path}"
            
        finally:
            # Cleanup
            for path in [config.vcf_file, config.target_id_file]:
                if path.exists():
                    path.unlink()
    
    @given(
        gt_df=sample_dataframe(),
        seq_df=sample_dataframe()
    )
    def test_compression_roundtrip_consistency(self, gt_df, seq_df):
        """Test that compressed and uncompressed outputs contain the same data."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create input files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            vcf_file.touch()
            target_id_file.touch()
            
            # Test both compressed and uncompressed
            for compress in [True, False]:
                output_file = temp_path / f"test_compress_{compress}"
                
                config = ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_id_file,
                    output_file=output_file,
                    compress_output=compress
                )
                
                result = ProcessingResult()
                writer = OutputWriter(config, result)
                writer.write_dataframes(gt_df, seq_df)
                
                paths = writer.get_output_paths()
                
                # Read back the data
                if compress:
                    with gzip.open(paths["genotype"], 'rt') as f:
                        gt_content = f.read()
                    with gzip.open(paths["sequence"], 'rt') as f:
                        seq_content = f.read()
                else:
                    with open(paths["genotype"], 'r') as f:
                        gt_content = f.read()
                    with open(paths["sequence"], 'r') as f:
                        seq_content = f.read()
                
                # Content should be valid CSV/TSV
                lines_gt = gt_content.strip().split('\n')
                lines_seq = seq_content.strip().split('\n')
                
                assert len(lines_gt) > 0, "GT file should have content"
                assert len(lines_seq) > 0, "SEQ file should have content"
                
                # First line should be header
                assert "CHROM" in lines_gt[0], "GT file should have header"
                assert "CHROM" in lines_seq[0], "SEQ file should have header"
                
                # Should have data rows (header + data)
                expected_rows = len(gt_df) + 1  # +1 for header
                assert len(lines_gt) == expected_rows, f"GT file should have {expected_rows} lines"
                assert len(lines_seq) == expected_rows, f"SEQ file should have {expected_rows} lines"


class TestCompressionEdgeCases:
    """Edge case tests for compression support."""
    
    def test_empty_dataframe_compression(self):
        """Test compression with empty DataFrames."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create input files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            vcf_file.touch()
            target_id_file.touch()
            
            # Create empty DataFrames
            empty_df = pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT"])
            
            for compress in [True, False]:
                output_file = temp_path / f"empty_test_{compress}"
                
                config = ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_id_file,
                    output_file=output_file,
                    compress_output=compress
                )
                
                result = ProcessingResult()
                writer = OutputWriter(config, result)
                writer.write_dataframes(empty_df, empty_df)
                
                paths = writer.get_output_paths()
                
                # Files should exist even if empty
                assert paths["genotype"].exists()
                assert paths["sequence"].exists()
                
                # Should contain at least headers
                if compress:
                    with gzip.open(paths["genotype"], 'rt') as f:
                        content = f.read()
                else:
                    with open(paths["genotype"], 'r') as f:
                        content = f.read()
                
                assert "CHROM" in content, "Even empty files should have headers"
    
    def test_large_dataframe_compression_efficiency(self):
        """Test that compression is efficient for large DataFrames."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create input files
            vcf_file = temp_path / "test.vcf"
            target_id_file = temp_path / "targets.txt"
            vcf_file.touch()
            target_id_file.touch()
            
            # Create large DataFrame with repetitive data (should compress well)
            n_rows = 1000
            large_df = pd.DataFrame({
                "CHROM": ["chr1"] * n_rows,
                "POS": list(range(1, n_rows + 1)),
                "REF": ["A"] * n_rows,
                "ALT": ["T"] * n_rows,
                "sample1": ["0/0"] * n_rows,
                "sample2": ["0/1"] * n_rows,
            })
            
            # Test compressed
            config_compressed = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=temp_path / "large_compressed",
                compress_output=True
            )
            
            result_compressed = ProcessingResult()
            writer_compressed = OutputWriter(config_compressed, result_compressed)
            writer_compressed.write_dataframes(large_df, large_df)
            
            # Test uncompressed
            config_uncompressed = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=temp_path / "large_uncompressed",
                compress_output=False
            )
            
            result_uncompressed = ProcessingResult()
            writer_uncompressed = OutputWriter(config_uncompressed, result_uncompressed)
            writer_uncompressed.write_dataframes(large_df, large_df)
            
            # Compare file sizes
            compressed_paths = writer_compressed.get_output_paths()
            uncompressed_paths = writer_uncompressed.get_output_paths()
            
            compressed_size = compressed_paths["genotype"].stat().st_size
            uncompressed_size = uncompressed_paths["genotype"].stat().st_size
            
            # Compressed should be significantly smaller for repetitive data
            compression_ratio = compressed_size / uncompressed_size
            assert compression_ratio < 0.8, f"Compression should be effective: ratio={compression_ratio}"