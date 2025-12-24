"""Property-based tests for configuration and data models."""

import tempfile
from pathlib import Path
from typing import List

import pytest
from hypothesis import given, strategies as st

from chip.vcf_processor.config import ProcessingConfig, ProcessingResult, VariantInfo


class TestVariantInfoProperties:
    """Property-based tests for VariantInfo dataclass."""
    
    @given(
        chrom=st.text(min_size=1, max_size=10).filter(lambda x: x.strip()),
        pos=st.integers(min_value=1, max_value=1000000),
        ref=st.text(alphabet="ATCG", min_size=1, max_size=10),
        alt=st.lists(st.text(alphabet="ATCG", min_size=1, max_size=10), min_size=1, max_size=5),
        genotypes=st.lists(st.sampled_from(["0/0", "0/1", "1/1", "./."]), min_size=1, max_size=10),
    )
    def test_variant_info_creation_and_validation(
        self, chrom: str, pos: int, ref: str, alt: List[str], genotypes: List[str]
    ) -> None:
        """Test that valid VariantInfo objects can be created and validated.
        
        **Feature: vcf-processor-optimization, Property 8: Parameter Validation**
        **Validates: Requirements 5.3**
        """
        variant_id = f"{chrom}_{pos}"
        
        # Valid variant should be created successfully
        variant = VariantInfo(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            variant_id=variant_id,
            genotypes=genotypes,
        )
        
        # Verify all fields are set correctly
        assert variant.chrom == chrom
        assert variant.pos == pos
        assert variant.ref == ref
        assert variant.alt == alt
        assert variant.variant_id == variant_id
        assert variant.genotypes == genotypes
        
        # Verify derived methods work
        assert variant.get_location_key() == variant_id
        
        variant_dict = variant.to_dict()
        assert variant_dict["CHROM"] == chrom
        assert variant_dict["POS"] == pos
        assert variant_dict["REF"] == ref
        assert variant_dict["ALT"] == ",".join(alt)
        assert variant_dict["ID"] == variant_id
    
    def test_variant_info_validation_errors(self) -> None:
        """Test that invalid VariantInfo parameters are rejected."""
        # Empty chromosome should raise ValueError
        with pytest.raises(ValueError, match="Chromosome cannot be empty"):
            VariantInfo(
                chrom="",
                pos=100,
                ref="A",
                alt=["T"],
                variant_id="test",
                genotypes=["0/1"],
            )
        
        # Zero or negative position should raise ValueError
        with pytest.raises(ValueError, match="Position must be positive"):
            VariantInfo(
                chrom="chr1",
                pos=0,
                ref="A",
                alt=["T"],
                variant_id="test",
                genotypes=["0/1"],
            )
        
        # Empty reference should raise ValueError
        with pytest.raises(ValueError, match="Reference allele cannot be empty"):
            VariantInfo(
                chrom="chr1",
                pos=100,
                ref="",
                alt=["T"],
                variant_id="test",
                genotypes=["0/1"],
            )
        
        # Empty alt list should raise ValueError
        with pytest.raises(ValueError, match="Alternative alleles cannot be empty"):
            VariantInfo(
                chrom="chr1",
                pos=100,
                ref="A",
                alt=[],
                variant_id="test",
                genotypes=["0/1"],
            )


class TestProcessingConfigProperties:
    """Property-based tests for ProcessingConfig dataclass."""
    
    @given(
        miss_fmt=st.text(min_size=1, max_size=5),
        gt_sep=st.text(max_size=2),
        threads=st.integers(min_value=1, max_value=16),
        batch_size=st.integers(min_value=1, max_value=100000),
        compress_output=st.booleans(),
        verbose=st.booleans(),
        quiet=st.booleans(),
        dry_run=st.booleans(),
    )
    def test_processing_config_validation(
        self,
        miss_fmt: str,
        gt_sep: str,
        threads: int,
        batch_size: int,
        compress_output: bool,
        verbose: bool,
        quiet: bool,
        dry_run: bool,
    ) -> None:
        """Test ProcessingConfig validation with various parameter combinations.
        
        **Feature: vcf-processor-optimization, Property 8: Parameter Validation**
        **Validates: Requirements 5.3**
        """
        # Skip invalid combinations
        if verbose and quiet:
            return
        
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            
            # Create dummy input files
            vcf_file = tmp_path / "test.vcf"
            vcf_file.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            
            target_file = tmp_path / "targets.txt"
            target_file.write_text("chr1_100\n")
            
            output_file = tmp_path / "output"
            
            # Valid configuration should be created successfully
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_file,
                output_file=output_file,
                miss_fmt=miss_fmt,
                gt_sep=gt_sep,
                threads=threads,
                batch_size=batch_size,
                compress_output=compress_output,
                verbose=verbose,
                quiet=quiet,
                dry_run=dry_run,
            )
            
            # Verify all fields are set correctly
            assert config.vcf_file == vcf_file
            assert config.target_id_file == target_file
            assert config.output_file == output_file
            assert config.miss_fmt == miss_fmt
            assert config.gt_sep == gt_sep
            assert config.threads == threads
            assert config.batch_size == batch_size
            assert config.compress_output == compress_output
            assert config.verbose == verbose
            assert config.quiet == quiet
            assert config.dry_run == dry_run
            
            # Verify to_dict method works
            config_dict = config.to_dict()
            assert config_dict["threads"] == threads
            assert config_dict["batch_size"] == batch_size
            assert config_dict["verbose"] == verbose
    
    def test_processing_config_validation_errors(self) -> None:
        """Test that invalid ProcessingConfig parameters are rejected."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            
            # Create valid files for most tests
            vcf_file = tmp_path / "test.vcf"
            vcf_file.write_text("##fileformat=VCFv4.2\n")
            
            target_file = tmp_path / "targets.txt"
            target_file.write_text("chr1_100\n")
            
            output_file = tmp_path / "output"
            
            # Non-existent VCF file should raise FileNotFoundError
            with pytest.raises(FileNotFoundError, match="VCF file not found"):
                ProcessingConfig(
                    vcf_file=tmp_path / "nonexistent.vcf",
                    target_id_file=target_file,
                    output_file=output_file,
                )
            
            # Non-existent target file should raise FileNotFoundError
            with pytest.raises(FileNotFoundError, match="Target ID file not found"):
                ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=tmp_path / "nonexistent.txt",
                    output_file=output_file,
                )
            
            # Zero threads should raise ValueError
            with pytest.raises(ValueError, match="Threads must be positive"):
                ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_file,
                    output_file=output_file,
                    threads=0,
                )
            
            # Zero batch size should raise ValueError
            with pytest.raises(ValueError, match="Batch size must be positive"):
                ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_file,
                    output_file=output_file,
                    batch_size=0,
                )
            
            # Empty miss_fmt should raise ValueError
            with pytest.raises(ValueError, match="Missing format cannot be empty"):
                ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_file,
                    output_file=output_file,
                    miss_fmt="",
                )
            
            # Both verbose and quiet should raise ValueError
            with pytest.raises(ValueError, match="Cannot enable both verbose and quiet modes"):
                ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_file,
                    output_file=output_file,
                    verbose=True,
                    quiet=True,
                )


class TestProcessingResultProperties:
    """Property-based tests for ProcessingResult dataclass."""
    
    @given(
        total_variants=st.integers(min_value=0, max_value=1000000),
        processed_variants=st.integers(min_value=0, max_value=1000000),
        skipped_variants=st.integers(min_value=0, max_value=1000000),
        processing_time=st.floats(min_value=0.0, max_value=3600.0),
    )
    def test_processing_result_calculations(
        self,
        total_variants: int,
        processed_variants: int,
        skipped_variants: int,
        processing_time: float,
    ) -> None:
        """Test ProcessingResult calculations and properties.
        
        **Feature: vcf-processor-optimization, Property 8: Parameter Validation**
        **Validates: Requirements 5.3**
        """
        # Ensure processed + skipped doesn't exceed total
        if processed_variants + skipped_variants > total_variants:
            processed_variants = min(processed_variants, total_variants)
            skipped_variants = total_variants - processed_variants
        
        result = ProcessingResult(
            total_variants=total_variants,
            processed_variants=processed_variants,
            skipped_variants=skipped_variants,
            processing_time=processing_time,
        )
        
        # Test success rate calculation
        if total_variants == 0:
            assert result.success_rate == 0.0
        else:
            expected_rate = processed_variants / total_variants
            assert abs(result.success_rate - expected_rate) < 1e-10
        
        # Test error handling
        assert not result.has_errors
        
        result.add_error("Test error")
        assert result.has_errors
        assert len(result.errors) == 1
        assert "Test error" in result.errors
        
        # Test summary generation
        summary = result.summary()
        assert "VCF Processing Summary" in summary
        assert f"{total_variants:,}" in summary  # Use comma formatting like the summary does
        assert f"{processed_variants:,}" in summary  # Use comma formatting for processed variants too
        
        # Test to_dict conversion
        result_dict = result.to_dict()
        assert result_dict["total_variants"] == total_variants
        assert result_dict["processed_variants"] == processed_variants
        assert result_dict["has_errors"] == result.has_errors