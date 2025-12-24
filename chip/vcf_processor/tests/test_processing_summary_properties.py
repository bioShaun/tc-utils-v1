"""Property-based tests for processing summary generation."""

import tempfile
from pathlib import Path
from typing import Any, Dict

import pytest
from hypothesis import given, strategies as st

from chip.vcf_processor.config import ProcessingConfig
from chip.vcf_processor.vcf_processor import VCFProcessor, ProcessingSummary


class TestProcessingSummaryProperties:
    """Property tests for processing summary generation."""
    
    @given(
        num_variants=st.integers(min_value=1, max_value=50),
        num_samples=st.integers(min_value=1, max_value=10),
        batch_size=st.integers(min_value=100, max_value=5000),
        threads=st.integers(min_value=1, max_value=4),
        compress_output=st.booleans(),
        verbose=st.booleans(),
        quiet=st.booleans(),
        dry_run=st.booleans()
    )
    def test_processing_summary_generation(
        self, 
        num_variants: int, 
        num_samples: int, 
        batch_size: int, 
        threads: int,
        compress_output: bool,
        verbose: bool,
        quiet: bool,
        dry_run: bool
    ):
        """Property 14: Processing Summary Generation
        
        Tests that processing summaries are generated correctly:
        - Summary contains all required fields
        - Summary data is consistent with processing results
        - Summary formatting is valid and readable
        - Summary reflects actual configuration and results
        """
        # Skip conflicting verbose/quiet combinations
        if verbose and quiet:
            return
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test VCF
            vcf_file = temp_path / "summary_test.vcf"
            vcf_lines = ["##fileformat=VCFv4.2"]
            vcf_lines.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
            
            # Create sample header
            sample_names = [f"sample{i+1}" for i in range(num_samples)]
            header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names)
            vcf_lines.append(header)
            
            # Generate variants
            target_ids = []
            chromosomes = ["chr1", "chr2", "chr3"]
            alleles = ["A", "T", "G", "C"]
            genotypes = ["0/0", "0/1", "1/1"]
            
            for i in range(num_variants):
                chrom = chromosomes[i % len(chromosomes)]
                pos = 1000 + i * 100
                ref = alleles[i % len(alleles)]
                alt = alleles[(i + 1) % len(alleles)]
                
                # Generate genotypes for all samples
                sample_gts = []
                for j in range(num_samples):
                    gt = genotypes[(i + j) % len(genotypes)]
                    sample_gts.append(gt)
                
                variant_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t.\tGT\t" + "\t".join(sample_gts)
                vcf_lines.append(variant_line)
                
                # Add to targets (include ~80% of variants)
                if i % 5 != 0:
                    target_ids.append(f"{chrom}_{pos}_{ref}_{alt}")
            
            vcf_file.write_text("\n".join(vcf_lines))
            
            # Create target file
            target_file = temp_path / "targets.txt"
            target_file.write_text("\n".join(target_ids))
            
            output_file = temp_path / "summary_output"
            
            # Create configuration
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_file,
                output_file=output_file,
                batch_size=batch_size,
                threads=threads,
                compress_output=compress_output,
                verbose=verbose,
                quiet=quiet,
                dry_run=dry_run
            )
            
            # Process VCF
            processor = VCFProcessor(config)
            result = processor.process()
            
            # Get processing summary
            summary_dict = processor.get_processing_summary()
            formatted_summary = processor.format_summary()
            
            # Verify summary contains required fields
            required_fields = [
                "input_vcf", "target_ids_file", "output_prefix",
                "variants_processed", "total_variants", "processing_rate",
                "batch_size", "threads", "compression_enabled", "dry_run",
                "elapsed_time_seconds", "output_files", "output_file_count",
                "errors", "error_count", "success", "processing_mode"
            ]
            
            for field in required_fields:
                assert field in summary_dict, f"Summary missing required field: {field}"
            
            # Verify summary data consistency
            assert summary_dict["input_vcf"] == str(config.vcf_file)
            assert summary_dict["target_ids_file"] == str(config.target_id_file)
            assert summary_dict["output_prefix"] == str(config.output_file)
            assert summary_dict["batch_size"] == config.batch_size
            assert summary_dict["threads"] == config.threads
            assert summary_dict["compression_enabled"] == config.compress_output
            assert summary_dict["dry_run"] == config.dry_run
            
            # Verify processing statistics consistency
            assert summary_dict["variants_processed"] == result.processed_variants
            assert summary_dict["variants_processed"] >= 0
            
            if summary_dict["total_variants"] > 0:
                assert summary_dict["variants_processed"] <= summary_dict["total_variants"]
            
            # Verify error information consistency
            assert summary_dict["error_count"] == len(result.errors)
            assert summary_dict["errors"] == result.errors
            assert summary_dict["success"] == (len(result.errors) == 0)
            
            # Verify output file information
            assert summary_dict["output_file_count"] == len(result.output_files)
            assert len(summary_dict["output_files"]) == len(result.output_files)
            
            # Verify timing information
            assert summary_dict["elapsed_time_seconds"] >= 0
            assert "start_time" in summary_dict
            
            if summary_dict["end_time"] is not None:
                assert summary_dict["end_time"] >= summary_dict["start_time"]
            
            # Verify processing rate calculation
            if summary_dict["elapsed_time_seconds"] > 0:
                expected_rate = summary_dict["variants_processed"] / summary_dict["elapsed_time_seconds"]
                assert abs(summary_dict["processing_rate"] - expected_rate) < 0.1
            else:
                assert summary_dict["processing_rate"] >= 0
            
            # Verify processing mode
            if summary_dict["batch_size"] >= 10000:
                assert summary_dict["processing_mode"] == "batch"
            else:
                assert summary_dict["processing_mode"] == "single_pass"
            
            # Verify formatted summary is valid
            assert isinstance(formatted_summary, str)
            assert len(formatted_summary) > 0
            
            # Verify formatted summary contains key information
            assert str(config.vcf_file) in formatted_summary
            assert str(summary_dict["variants_processed"]) in formatted_summary
            
            if summary_dict["success"]:
                assert "✓" in formatted_summary or "completed successfully" in formatted_summary.lower()
            else:
                assert "✗" in formatted_summary or "error" in formatted_summary.lower()
            
            # Verify summary sections are present
            expected_sections = ["INPUT/OUTPUT", "PROCESSING STATISTICS", "CONFIGURATION", "TIMING", "STATUS"]
            for section in expected_sections:
                assert section in formatted_summary, f"Formatted summary missing section: {section}"
    
    @given(
        processing_errors=st.lists(st.text(min_size=1, max_size=100), min_size=0, max_size=5),
        variants_processed=st.integers(min_value=0, max_value=1000),
        total_variants=st.integers(min_value=0, max_value=1000)
    )
    def test_summary_error_handling(
        self, 
        processing_errors: list[str], 
        variants_processed: int, 
        total_variants: int
    ):
        """Test summary generation with various error conditions."""
        # Ensure total_variants >= variants_processed
        if total_variants < variants_processed:
            total_variants = variants_processed
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create minimal test files
            vcf_file = temp_path / "error_test.vcf"
            vcf_file.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n")
            
            target_file = temp_path / "targets.txt"
            target_file.write_text("dummy_target\n")
            
            output_file = temp_path / "error_output"
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_file,
                output_file=output_file,
                dry_run=True,  # Use dry run to avoid complex processing
                quiet=True
            )
            
            # Create processor and manually set result data
            processor = VCFProcessor(config)
            
            # Manually set processing results to test summary generation
            processor.result.processed_variants = variants_processed
            processor.result.total_variants = total_variants
            
            # Add errors to result
            for error in processing_errors:
                processor.result.add_error(error)
            
            # Finalize summary
            processor.summary.finalize()
            
            # Get summary
            summary_dict = processor.get_processing_summary()
            formatted_summary = processor.format_summary()
            
            # Verify error handling in summary
            assert summary_dict["error_count"] == len(processing_errors)
            assert summary_dict["errors"] == processing_errors
            assert summary_dict["success"] == (len(processing_errors) == 0)
            
            # Verify variants information
            assert summary_dict["variants_processed"] == variants_processed
            assert summary_dict["total_variants"] == total_variants
            
            # Verify error information in formatted summary
            if processing_errors:
                assert "ERRORS:" in formatted_summary
                assert "✗" in formatted_summary or "error" in formatted_summary.lower()
                
                # Should show first few errors
                for error in processing_errors[:3]:  # First 3 errors should be shown
                    assert error in formatted_summary
            else:
                assert "✓" in formatted_summary or "successfully" in formatted_summary.lower()
    
    @given(
        file_sizes=st.lists(st.integers(min_value=1024, max_value=10*1024*1024), min_size=0, max_size=3),
        elapsed_time=st.floats(min_value=0.1, max_value=3600.0),
        processing_mode=st.sampled_from(["batch", "single_pass"])
    )
    def test_summary_performance_metrics(
        self, 
        file_sizes: list[int], 
        elapsed_time: float, 
        processing_mode: str
    ):
        """Test summary generation with performance metrics."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test files
            vcf_file = temp_path / "perf_test.vcf"
            vcf_file.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n")
            
            target_file = temp_path / "targets.txt"
            target_file.write_text("dummy_target\n")
            
            # Create output files with specified sizes (simulate processing output)
            output_files = []
            for i, size in enumerate(file_sizes):
                output_file = temp_path / f"output_{i}.txt"
                output_file.write_bytes(b"x" * size)  # Create file with specified size
                output_files.append(output_file)
            
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_file,
                output_file=temp_path / "perf_output",
                batch_size=10000 if processing_mode == "batch" else 1000,
                dry_run=True,
                quiet=True
            )
            
            # Create processor
            processor = VCFProcessor(config)
            
            # Manually set timing and output files
            processor.summary.start_time = 1000.0  # Fixed start time
            processor.summary.end_time = processor.summary.start_time + elapsed_time
            
            # Add output files to result
            for output_file in output_files:
                processor.result.output_files.append(output_file)
            
            # Set some processing statistics
            variants_processed = int(elapsed_time * 10)  # Simulate processing rate
            processor.result.processed_variants = variants_processed
            processor.result.total_variants = variants_processed
            
            # Get summary
            summary_dict = processor.get_processing_summary()
            formatted_summary = processor.format_summary()
            
            # Verify performance metrics
            assert abs(summary_dict["elapsed_time_seconds"] - elapsed_time) < 0.1
            
            expected_rate = variants_processed / elapsed_time
            assert abs(summary_dict["processing_rate"] - expected_rate) < 0.1
            
            # Verify file size information
            if file_sizes:
                assert "total_output_size_mb" in summary_dict
                expected_total_size = sum(file_sizes) / (1024 * 1024)
                assert abs(summary_dict["total_output_size_mb"] - expected_total_size) < 0.1
            
            # Verify processing mode
            assert summary_dict["processing_mode"] == processing_mode
            
            # Verify formatted summary includes performance info
            assert f"{elapsed_time:.2f} seconds" in formatted_summary
            assert f"{expected_rate:.1f} variants/sec" in formatted_summary
            
            if elapsed_time >= 60:
                # Should show minutes and seconds for long runs
                minutes = int(elapsed_time // 60)
                assert f"{minutes}m" in formatted_summary
            
            # Verify output file information
            assert summary_dict["output_file_count"] == len(output_files)
            
            if output_files:
                assert "OUTPUT FILES:" in formatted_summary
                for output_file in output_files:
                    assert str(output_file) in formatted_summary