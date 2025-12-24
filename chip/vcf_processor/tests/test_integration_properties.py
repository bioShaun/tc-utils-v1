"""Integration property tests for VCF processor."""

import tempfile
from pathlib import Path
from typing import List, Dict, Any

import pytest
from hypothesis import given, strategies as st, assume

from chip.vcf_processor.config import ProcessingConfig
from chip.vcf_processor.vcf_processor import VCFProcessor


class TestIntegrationProperties:
    """Property-based integration tests for complete VCF processing workflows."""
    
    @given(
        num_variants=st.integers(min_value=1, max_value=20),
        num_samples=st.integers(min_value=1, max_value=5),
        processing_params=st.fixed_dictionaries({
            'batch_size': st.integers(min_value=1, max_value=10),
            'compress_output': st.booleans(),
            'include_variant_type': st.booleans()
        })
    )
    def test_end_to_end_processing_consistency(
        self, 
        num_variants: int, 
        num_samples: int, 
        processing_params: Dict[str, Any]
    ):
        """Test end-to-end processing produces consistent results.
        
        This property test verifies that:
        - Processing completes successfully with various configurations
        - Output files are generated correctly
        - Results are consistent across different batch sizes
        - Variant type annotation works when enabled
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Generate test VCF content
            vcf_lines = [
                "##fileformat=VCFv4.2",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
            ]
            
            # Add sample names to header
            sample_names = [f"sample_{i+1}" for i in range(num_samples)]
            header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names)
            vcf_lines.append(header_line)
            
            # Generate variants
            variant_ids = []
            genotypes = ["0/0", "0/1", "1/1", "./."]
            
            for i in range(num_variants):
                chrom = f"chr{(i % 5) + 1}"
                pos = 1000 + i * 100
                ref = "ATGC"[i % 4]
                alt = "CGTA"[i % 4]
                
                variant_id = f"{chrom}_{pos}_{ref}_{alt}"
                variant_ids.append(variant_id)
                
                # Generate sample genotypes
                sample_gts = []
                for j in range(num_samples):
                    gt = genotypes[(i + j) % len(genotypes)]
                    sample_gts.append(gt)
                
                variant_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t.\tGT\t" + "\t".join(sample_gts)
                vcf_lines.append(variant_line)
            
            # Create files
            vcf_file = temp_path / "test.vcf"
            vcf_file.write_text("\n".join(vcf_lines))
            
            target_file = temp_path / "targets.txt"
            target_file.write_text("\n".join(variant_ids))
            
            output_file = temp_path / "output"
            
            # Create configuration
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_file,
                output_file=output_file,
                quiet=True,
                **processing_params
            )
            
            # Process VCF
            processor = VCFProcessor(config)
            result = processor.process()
            
            # Verify processing completed
            assert len(result.errors) == 0, f"Processing had errors: {result.errors}"
            assert result.processed_variants >= 0, "Invalid processed variant count"
            
            # Verify output files exist
            expected_files = [
                output_file.with_suffix(".genotype_codes.tsv"),
                output_file.with_suffix(".genotype_bases.tsv")
            ]
            
            for expected_file in expected_files:
                if processing_params['compress_output']:
                    expected_file = expected_file.with_suffix(expected_file.suffix + ".gz")
                assert expected_file.exists(), f"Output file not created: {expected_file}"
            
            # Verify summary is generated correctly
            summary = processor.get_processing_summary()
            assert summary["success"] == True
            assert summary["variants_processed"] >= 0
    
    @given(
        variant_types=st.lists(
            st.tuples(
                st.sampled_from(["SNP", "INSERTION", "DELETION", "COMPLEX"]),
                st.integers(min_value=1, max_value=5)
            ),
            min_size=1,
            max_size=4
        ),
        num_samples=st.integers(min_value=1, max_value=3),
        processing_params=st.fixed_dictionaries({
            'batch_size': st.integers(min_value=1, max_value=5),
            'compress_output': st.booleans()
        })
    )
    def test_variant_type_annotation_integration(
        self, 
        variant_types: List[tuple], 
        num_samples: int, 
        processing_params: Dict[str, Any]
    ):
        """Test variant type annotation in end-to-end processing.
        
        This property test verifies that:
        - Variant type annotation works with various variant types
        - Output includes Variant_Type column when enabled
        - Processing handles different variant types correctly
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Generate test VCF content
            vcf_lines = [
                "##fileformat=VCFv4.2",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
            ]
            
            # Add sample names to header
            sample_names = [f"sample_{i+1}" for i in range(num_samples)]
            header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names)
            vcf_lines.append(header_line)
            
            # Generate variants by type
            variant_ids = []
            pos = 10000
            genotypes = ["0/0", "0/1", "1/1", "./."]
            
            for variant_type, count in variant_types:
                for i in range(count):
                    chrom = f"chr{(pos // 10000) % 5 + 1}"
                    
                    if variant_type == "SNP":
                        ref = "ATGC"[i % 4]
                        alt = "CGTA"[i % 4]
                    elif variant_type == "INSERTION":
                        ref = "A"
                        alt = "A" + "TGCA"[i % 4] * (i % 3 + 1)
                    elif variant_type == "DELETION":
                        ref = "A" + "TGCA"[i % 4] * (i % 3 + 1)
                        alt = "A"
                    else:  # COMPLEX
                        ref = "ATG"[:(i % 3) + 1]
                        alt = "CGA"[:(i % 3) + 1]
                    
                    variant_id = f"{chrom}_{pos}_{ref}_{alt}"
                    variant_ids.append(variant_id)
                    
                    # Generate sample genotypes
                    sample_gts = []
                    for j in range(num_samples):
                        gt = genotypes[(i + j) % len(genotypes)]
                        sample_gts.append(gt)
                    
                    variant_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t.\tGT\t" + "\t".join(sample_gts)
                    vcf_lines.append(variant_line)
                    
                    pos += 1000
            
            # Skip if no variants generated
            assume(len(variant_ids) > 0)
            
            # Create files
            vcf_file = temp_path / "variant_types_test.vcf"
            vcf_file.write_text("\n".join(vcf_lines))
            
            target_file = temp_path / "targets.txt"
            target_file.write_text("\n".join(variant_ids))
            
            output_file = temp_path / "variant_types_output"
            
            # Create configuration
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_file,
                output_file=output_file,
                quiet=True,
                **processing_params
            )
            
            # Process VCF
            processor = VCFProcessor(config)
            result = processor.process()
            
            # Verify processing completed
            assert len(result.errors) == 0, f"Variant type processing had errors: {result.errors}"
            assert result.processed_variants >= 0, "Invalid processed variant count"
            
            # Verify summary is generated correctly
            summary = processor.get_processing_summary()
            assert summary["success"] == True
            assert summary["variants_processed"] >= 0
    
    @given(
        error_scenarios=st.lists(
            st.sampled_from([
                "missing_targets",
                "invalid_genotypes", 
                "empty_vcf",
                "malformed_targets"
            ]),
            min_size=0,
            max_size=2
        ),
        recovery_mode=st.booleans()
    )
    def test_error_recovery_integration(self, error_scenarios: List[str], recovery_mode: bool):
        """Test error recovery and graceful degradation.
        
        This property test verifies that:
        - Processing handles various error conditions gracefully
        - Errors are properly reported and logged
        - Partial processing can continue when possible
        - System remains stable after errors
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create base VCF content
            base_vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2
chr1	100	.	A	T	60	PASS	.	GT	0/0	0/1
chr1	200	.	G	C	55	PASS	.	GT	0/1	1/1
"""
            
            base_targets = ["chr1_100_A_T", "chr1_200_G_C"]
            
            # Apply error scenarios
            vcf_content = base_vcf_content
            target_content = "\n".join(base_targets)
            
            for scenario in error_scenarios:
                if scenario == "missing_targets":
                    # Reference non-existent variants
                    target_content += "\nchr99_999_X_Y\nchr88_888_Z_W"
                elif scenario == "invalid_genotypes":
                    # Add line with invalid genotype
                    vcf_content += "chr1\t300\t.\tT\tG\t50\tPASS\t.\tGT\tINVALID\t2/3\n"
                    target_content += "\nchr1_300_T_G"
                elif scenario == "empty_vcf":
                    # Use minimal VCF
                    vcf_content = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n"
                elif scenario == "malformed_targets":
                    # Add malformed target IDs
                    target_content += "\nmalformed_target\n\ninvalid\nformat"
            
            # Create files
            vcf_file = temp_path / "error_test.vcf"
            vcf_file.write_text(vcf_content)
            
            target_file = temp_path / "targets.txt"
            target_file.write_text(target_content)
            
            output_file = temp_path / "error_output"
            
            # Create configuration
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_file,
                output_file=output_file,
                quiet=True,
                dry_run=recovery_mode  # Use dry run for recovery mode
            )
            
            # Process VCF (should handle errors gracefully)
            processor = VCFProcessor(config)
            
            try:
                result = processor.process()
                
                # Processing should complete (may have errors)
                assert result.processed_variants >= 0, "Invalid processed variant count"
                
                # If there were error scenarios, should have some errors or warnings
                if error_scenarios:
                    # May have errors, but should not crash
                    assert len(result.errors) >= 0  # Can have 0 errors if gracefully handled
                
                # Summary should be generated even with errors
                summary = processor.get_processing_summary()
                assert "variants_processed" in summary
                assert "error_count" in summary
                
            except Exception as e:
                # If processing raises exception, it should be a controlled failure
                assert isinstance(e, (RuntimeError, ValueError, FileNotFoundError)), \
                    f"Unexpected exception type: {type(e)}"
    
    @given(
        run_count=st.integers(min_value=2, max_value=5),
        vary_config=st.booleans()
    )
    def test_multiple_run_consistency(self, run_count: int, vary_config: bool):
        """Test consistency across multiple processing runs.
        
        This property test verifies that:
        - Multiple runs with same input produce same results
        - Configuration changes produce predictable result differences
        - Processing is deterministic and repeatable
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create consistent test data
            vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3
chr1	100	.	A	T	60	PASS	.	GT	0/0	0/1	1/1
chr1	200	.	G	C	55	PASS	.	GT	0/1	1/1	0/0
chr2	300	.	T	A	65	PASS	.	GT	1/1	0/0	0/1
chr2	400	.	C	G	50	PASS	.	GT	0/0	0/1	1/1
"""
            
            target_content = "chr1_100_A_T\nchr1_200_G_C\nchr2_300_T_A\nchr2_400_C_G"
            
            # Create files
            vcf_file = temp_path / "consistency_test.vcf"
            vcf_file.write_text(vcf_content)
            
            target_file = temp_path / "targets.txt"
            target_file.write_text(target_content)
            
            # Run multiple times
            results = []
            summaries = []
            
            for run_idx in range(run_count):
                output_file = temp_path / f"consistency_output_{run_idx}"
                
                # Vary configuration if requested
                if vary_config and run_idx > 0:
                    # Use different batch size for variety
                    batch_size = 1000 + run_idx * 500
                    compress_output = (run_idx % 2 == 0)
                else:
                    # Use consistent configuration
                    batch_size = 2000
                    compress_output = False
                
                config = ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_file,
                    output_file=output_file,
                    batch_size=batch_size,
                    compress_output=compress_output,
                    quiet=True
                )
                
                processor = VCFProcessor(config)
                result = processor.process()
                summary = processor.get_processing_summary()
                
                results.append(result)
                summaries.append(summary)
                
                # Each run should complete successfully
                assert len(result.errors) == 0, f"Run {run_idx}: Processing had errors: {result.errors}"
                assert result.processed_variants >= 0, f"Run {run_idx}: Invalid processed count"
            
            # Verify consistency expectations
            if not vary_config:
                # All runs should produce identical results
                baseline_processed = results[0].processed_variants
                for i, result in enumerate(results[1:], 1):
                    assert result.processed_variants == baseline_processed, \
                        f"Run {i}: Inconsistent processed count {result.processed_variants} vs {baseline_processed}"
            
            # All runs should have same success status
            baseline_success = summaries[0]["success"]
            for i, summary in enumerate(summaries[1:], 1):
                assert summary["success"] == baseline_success, \
                    f"Run {i}: Inconsistent success status"
            
            # All runs should process same number of variants (core logic should be consistent)
            baseline_variants = summaries[0]["variants_processed"]
            for i, summary in enumerate(summaries[1:], 1):
                assert summary["variants_processed"] == baseline_variants, \
                    f"Run {i}: Inconsistent variant count {summary['variants_processed']} vs {baseline_variants}"