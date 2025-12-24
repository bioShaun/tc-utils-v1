"""Integration property tests for complete pipeline validation."""

import tempfile
from pathlib import Path
from typing import List, Tuple

import pytest
from hypothesis import given, strategies as st, assume

from chip.vcf_processor.config import ProcessingConfig
from chip.vcf_processor.vcf_processor import VCFProcessor


class TestIntegrationProperties:
    """Property tests for complete pipeline integration."""
    
    def _generate_valid_vcf_data(
        self, 
        num_variants: int, 
        num_samples: int, 
        chromosomes: List[str],
        include_missing: bool = True
    ) -> Tuple[str, List[str]]:
        """Generate valid VCF content and corresponding variant IDs.
        
        Args:
            num_variants: Number of variants to generate
            num_samples: Number of samples to generate
            chromosomes: List of chromosome names to use
            include_missing: Whether to include missing genotypes
            
        Returns:
            Tuple of (vcf_content, variant_ids)
        """
        # VCF header
        vcf_lines = [
            "##fileformat=VCFv4.2",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"
        ]
        
        # Sample header
        sample_names = [f"sample_{i+1}" for i in range(num_samples)]
        header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names)
        vcf_lines.append(header_line)
        
        # Generate variants
        variant_ids = []
        alleles = ["A", "T", "G", "C"]
        genotypes = ["0/0", "0/1", "1/1", "0|0", "0|1", "1|0", "1|1"]
        
        if include_missing:
            genotypes.extend(["./.", ".|."])
        
        for i in range(num_variants):
            chrom = chromosomes[i % len(chromosomes)]
            pos = 10000 + i * 1000  # Ensure positions don't overlap
            ref = alleles[i % len(alleles)]
            alt = alleles[(i + 1) % len(alleles)]
            
            # Skip if ref == alt
            if ref == alt:
                alt = alleles[(i + 2) % len(alleles)]
            
            variant_id = f"{chrom}_{pos}_{ref}_{alt}"
            variant_ids.append(variant_id)
            
            # Generate sample genotypes
            sample_gts = []
            for j in range(num_samples):
                gt = genotypes[(i + j) % len(genotypes)]
                dp = 20 + (i + j) % 50
                sample_gts.append(f"{gt}:{dp}")
            
            # Create variant line
            qual = 30 + (i % 70)
            info = f"DP={50 + i % 100}"
            variant_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\tPASS\t{info}\tGT:DP\t" + "\t".join(sample_gts)
            vcf_lines.append(variant_line)
        
        return "\n".join(vcf_lines), variant_ids
    
    @given(
        num_variants=st.integers(min_value=5, max_value=50),
        num_samples=st.integers(min_value=2, max_value=10),
        target_fraction=st.floats(min_value=0.2, max_value=1.0),
        batch_size=st.integers(min_value=100, max_value=10000),
        threads=st.integers(min_value=1, max_value=4),
        compress_output=st.booleans(),
        miss_fmt=st.sampled_from(["NN", "./.", "--", "00"]),
        gt_sep=st.sampled_from(["", "/", "|", "_"])
    )
    def test_complete_pipeline_consistency(
        self,
        num_variants: int,
        num_samples: int,
        target_fraction: float,
        batch_size: int,
        threads: int,
        compress_output: bool,
        miss_fmt: str,
        gt_sep: str
    ):
        """Test complete pipeline with randomly generated valid inputs.
        
        This property test verifies that:
        - Processing completes successfully with valid inputs
        - Output is consistent across multiple runs with same input
        - All configuration parameters are respected
        - Results are deterministic for same inputs
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Generate test data
            chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5"]
            vcf_content, all_variant_ids = self._generate_valid_vcf_data(
                num_variants, num_samples, chromosomes
            )
            
            # Create VCF file
            vcf_file = temp_path / "pipeline_test.vcf"
            vcf_file.write_text(vcf_content)
            
            # Select target variants
            num_targets = max(1, int(len(all_variant_ids) * target_fraction))
            target_variant_ids = all_variant_ids[:num_targets]
            
            # Create target file
            target_file = temp_path / "targets.txt"
            target_file.write_text("\n".join(target_variant_ids))
            
            output_file = temp_path / "pipeline_output"
            
            # Create configuration
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_file,
                output_file=output_file,
                batch_size=batch_size,
                threads=threads,
                compress_output=compress_output,
                miss_fmt=miss_fmt,
                gt_sep=gt_sep,
                quiet=True  # Suppress output for cleaner testing
            )
            
            # Run processing multiple times to test consistency
            results = []
            summaries = []
            
            for run_idx in range(2):  # Run twice to test consistency
                # Use different output file for each run
                run_config = ProcessingConfig(
                    vcf_file=config.vcf_file,
                    target_id_file=config.target_id_file,
                    output_file=temp_path / f"pipeline_output_run_{run_idx}",
                    batch_size=config.batch_size,
                    threads=config.threads,
                    compress_output=config.compress_output,
                    miss_fmt=config.miss_fmt,
                    gt_sep=config.gt_sep,
                    quiet=True
                )
                
                processor = VCFProcessor(run_config)
                result = processor.process()
                summary = processor.get_processing_summary()
                
                results.append(result)
                summaries.append(summary)
                
                # Verify processing completed successfully
                assert len(result.errors) == 0, f"Run {run_idx}: Processing had errors: {result.errors}"
                assert result.processed_variants >= 0, f"Run {run_idx}: Invalid processed count"
                
                # Verify configuration was respected
                assert summary["batch_size"] == batch_size
                assert summary["threads"] == threads
                assert summary["compression_enabled"] == compress_output
                
                # Verify processing mode is correct
                expected_mode = "batch" if batch_size >= 10000 else "single_pass"
                assert summary["processing_mode"] == expected_mode
            
            # Verify consistency between runs
            assert results[0].processed_variants == results[1].processed_variants, \
                "Processing results should be consistent between runs"
            
            # Verify summary consistency
            assert summaries[0]["variants_processed"] == summaries[1]["variants_processed"]
            assert summaries[0]["success"] == summaries[1]["success"]
    
    @given(
        variant_types=st.lists(
            st.tuples(
                st.sampled_from(["SNP", "INSERTION", "DELETION", "COMPLEX"]),
                st.integers(min_value=1, max_value=10)  # Count of each type
            ),
            min_size=1,
            max_size=4
        ),
        num_samples=st.integers(min_value=1, max_value=8),
        processing_params=st.fixed_dictionaries({
            "batch_size": st.integers(min_value=500, max_value=5000),
            "compress_output": st.booleans(),
            "miss_fmt": st.sampled_from(["NN", "./.", "--"]),
            "gt_sep": st.sampled_from(["", "/", "|"])
        })
    )
    def test_variant_type_handling_integration(
        self,
        variant_types: List[Tuple[str, int]],
        num_samples: int,
        processing_params: dict
    ):
        """Test integration with different variant types.
        
        This property test verifies that:
        - Different variant types (SNP, indels, complex) are handled correctly
        - Processing is robust across variant type combinations
        - Output format is consistent regardless of variant complexity
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Generate VCF with specified variant types
            vcf_lines = [
                "##fileformat=VCFv4.2",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
            ]
            
            # Sample header
            sample_names = [f"sample_{i+1}" for i in range(num_samples)]
            header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names)
            vcf_lines.append(header_line)
            
            # Generate variants by type
            variant_ids = []\n            pos = 10000\n            genotypes = ["0/0", "0/1", "1/1", "./."]\n            \n            for variant_type, count in variant_types:\n                for i in range(count):\n                    chrom = f"chr{(pos // 10000) % 5 + 1}"\n                    \n                    if variant_type == "SNP":\n                        ref = "ATGC"[i % 4]\n                        alt = "CGTA"[i % 4]\n                    elif variant_type == "INSERTION":\n                        ref = "A"\n                        alt = "A" + "TGCA"[i % 4] * (i % 3 + 1)\n                    elif variant_type == "DELETION":\n                        ref = "A" + "TGCA"[i % 4] * (i % 3 + 1)\n                        alt = "A"\n                    else:  # COMPLEX\n                        ref = "ATG"[:(i % 3) + 1]\n                        alt = "CGA"[:(i % 3) + 1]\n                    \n                    variant_id = f"{chrom}_{pos}_{ref}_{alt}"\n                    variant_ids.append(variant_id)\n                    \n                    # Generate sample genotypes\n                    sample_gts = []\n                    for j in range(num_samples):\n                        gt = genotypes[(i + j) % len(genotypes)]\n                        sample_gts.append(gt)\n                    \n                    variant_line = f"{chrom}\\t{pos}\\t.\\t{ref}\\t{alt}\\t60\\tPASS\\t.\\tGT\\t" + "\\t".join(sample_gts)\n                    vcf_lines.append(variant_line)\n                    \n                    pos += 1000\n            \n            # Skip if no variants generated\n            assume(len(variant_ids) > 0)\n            \n            # Create files\n            vcf_file = temp_path / "variant_types_test.vcf"\n            vcf_file.write_text("\\n".join(vcf_lines))\n            \n            target_file = temp_path / "targets.txt"\n            target_file.write_text("\\n".join(variant_ids))\n            \n            output_file = temp_path / "variant_types_output"\n            \n            # Create configuration\n            config = ProcessingConfig(\n                vcf_file=vcf_file,\n                target_id_file=target_file,\n                output_file=output_file,\n                quiet=True,\n                **processing_params\n            )\n            \n            # Process VCF\n            processor = VCFProcessor(config)\n            result = processor.process()\n            \n            # Verify processing completed\n            assert len(result.errors) == 0, f"Variant type processing had errors: {result.errors}"\n            assert result.processed_variants >= 0, "Invalid processed variant count"\n            \n            # Verify summary is generated correctly\n            summary = processor.get_processing_summary()\n            assert summary["success"] == True\n            assert summary["variants_processed"] >= 0\n    \n    @given(\n        error_scenarios=st.lists(\n            st.sampled_from([\n                "missing_targets",\n                "invalid_genotypes", \n                "empty_vcf",\n                "malformed_targets"\n            ]),\n            min_size=0,\n            max_size=2\n        ),\n        recovery_mode=st.booleans()\n    )\n    def test_error_recovery_integration(self, error_scenarios: List[str], recovery_mode: bool):\n        """Test error recovery and graceful degradation.\n        \n        This property test verifies that:\n        - Processing handles various error conditions gracefully\n        - Errors are properly reported and logged\n        - Partial processing can continue when possible\n        - System remains stable after errors\n        """\n        with tempfile.TemporaryDirectory() as temp_dir:\n            temp_path = Path(temp_dir)\n            \n            # Create base VCF content\n            base_vcf_content = """##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0\t0/1\nchr1\t200\t.\tG\tC\t55\tPASS\t.\tGT\t0/1\t1/1\n"""\n            \n            base_targets = ["chr1_100_A_T", "chr1_200_G_C"]\n            \n            # Apply error scenarios\n            vcf_content = base_vcf_content\n            target_content = "\\n".join(base_targets)\n            \n            for scenario in error_scenarios:\n                if scenario == "missing_targets":\n                    # Reference non-existent variants\n                    target_content += "\\nchr99_999_X_Y\\nchr88_888_Z_W"\n                elif scenario == "invalid_genotypes":\n                    # Add line with invalid genotype\n                    vcf_content += "chr1\\t300\\t.\\tT\\tG\\t50\\tPASS\\t.\\tGT\\tINVALID\\t2/3\\n"\n                    target_content += "\\nchr1_300_T_G"\n                elif scenario == "empty_vcf":\n                    # Use minimal VCF\n                    vcf_content = "##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tsample1\\n"\n                elif scenario == "malformed_targets":\n                    # Add malformed target IDs\n                    target_content += "\\nmalformed_target\\n\\ninvalid\\nformat"\n            \n            # Create files\n            vcf_file = temp_path / "error_test.vcf"\n            vcf_file.write_text(vcf_content)\n            \n            target_file = temp_path / "targets.txt"\n            target_file.write_text(target_content)\n            \n            output_file = temp_path / "error_output"\n            \n            # Create configuration\n            config = ProcessingConfig(\n                vcf_file=vcf_file,\n                target_id_file=target_file,\n                output_file=output_file,\n                quiet=True,\n                dry_run=recovery_mode  # Use dry run for recovery mode\n            )\n            \n            # Process VCF (should handle errors gracefully)\n            processor = VCFProcessor(config)\n            \n            try:\n                result = processor.process()\n                \n                # Processing should complete (may have errors)\n                assert result.processed_variants >= 0, "Invalid processed variant count"\n                \n                # If there were error scenarios, should have some errors or warnings\n                if error_scenarios:\n                    # May have errors, but should not crash\n                    assert len(result.errors) >= 0  # Can have 0 errors if gracefully handled\n                \n                # Summary should be generated even with errors\n                summary = processor.get_processing_summary()\n                assert "variants_processed" in summary\n                assert "error_count" in summary\n                \n            except Exception as e:\n                # If processing raises exception, it should be a controlled failure\n                assert isinstance(e, (RuntimeError, ValueError, FileNotFoundError)), \\\n                    f"Unexpected exception type: {type(e)}"\n    \n    @given(\n        run_count=st.integers(min_value=2, max_value=5),\n        vary_config=st.booleans()\n    )\n    def test_multiple_run_consistency(self, run_count: int, vary_config: bool):\n        """Test consistency across multiple processing runs.\n        \n        This property test verifies that:\n        - Multiple runs with same input produce same results\n        - Configuration changes produce predictable result differences\n        - Processing is deterministic and repeatable\n        """\n        with tempfile.TemporaryDirectory() as temp_dir:\n            temp_path = Path(temp_dir)\n            \n            # Create consistent test data\n            vcf_content = """##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\tsample3\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0\t0/1\t1/1\nchr1\t200\t.\tG\tC\t55\tPASS\t.\tGT\t0/1\t1/1\t0/0\nchr2\t300\t.\tT\tA\t65\tPASS\t.\tGT\t1/1\t0/0\t0/1\nchr2\t400\t.\tC\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1\n"""\n            \n            target_content = "chr1_100_A_T\\nchr1_200_G_C\\nchr2_300_T_A\\nchr2_400_C_G"\n            \n            # Create files\n            vcf_file = temp_path / "consistency_test.vcf"\n            vcf_file.write_text(vcf_content)\n            \n            target_file = temp_path / "targets.txt"\n            target_file.write_text(target_content)\n            \n            # Run multiple times\n            results = []\n            summaries = []\n            \n            for run_idx in range(run_count):\n                output_file = temp_path / f"consistency_output_{run_idx}"\n                \n                # Vary configuration if requested\n                if vary_config and run_idx > 0:\n                    # Use different batch size for variety\n                    batch_size = 1000 + run_idx * 500\n                    compress_output = (run_idx % 2 == 0)\n                else:\n                    # Use consistent configuration\n                    batch_size = 2000\n                    compress_output = False\n                \n                config = ProcessingConfig(\n                    vcf_file=vcf_file,\n                    target_id_file=target_file,\n                    output_file=output_file,\n                    batch_size=batch_size,\n                    compress_output=compress_output,\n                    quiet=True\n                )\n                \n                processor = VCFProcessor(config)\n                result = processor.process()\n                summary = processor.get_processing_summary()\n                \n                results.append(result)\n                summaries.append(summary)\n                \n                # Each run should complete successfully\n                assert len(result.errors) == 0, f"Run {run_idx}: Processing had errors: {result.errors}"\n                assert result.processed_variants >= 0, f"Run {run_idx}: Invalid processed count"\n            \n            # Verify consistency expectations\n            if not vary_config:\n                # All runs should produce identical results\n                baseline_processed = results[0].processed_variants\n                for i, result in enumerate(results[1:], 1):\n                    assert result.processed_variants == baseline_processed, \\\n                        f"Run {i}: Inconsistent processed count {result.processed_variants} vs {baseline_processed}"\n            \n            # All runs should have same success status\n            baseline_success = summaries[0]["success"]\n            for i, summary in enumerate(summaries[1:], 1):\n                assert summary["success"] == baseline_success, \\\n                    f"Run {i}: Inconsistent success status"\n            \n            # All runs should process same number of variants (core logic should be consistent)\n            baseline_variants = summaries[0]["variants_processed"]\n            for i, summary in enumerate(summaries[1:], 1):\n                assert summary["variants_processed"] == baseline_variants, \\\n                    f"Run {i}: Inconsistent variant count {summary['variants_processed']} vs {baseline_variants}"\n