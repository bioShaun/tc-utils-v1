"""Configuration and data models for VCF processor."""

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

from loguru import logger


@dataclass
class VariantInfo:
    """Information about a VCF variant.
    
    Attributes:
        chrom: Chromosome name
        pos: Position (1-based)
        ref: Reference allele
        alt: List of alternative alleles
        variant_id: Unique variant identifier (CHROM_POS format)
        genotypes: List of raw GT values for all samples
    """
    chrom: str
    pos: int
    ref: str
    alt: List[str]
    variant_id: str
    genotypes: List[str]
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            "CHROM": self.chrom,
            "POS": self.pos,
            "REF": self.ref,
            "ALT": ",".join(self.alt) if self.alt else ".",
            "ID": self.variant_id,
            "genotypes": self.genotypes,
        }
    
    def get_location_key(self) -> str:
        """Get location key for variant."""
        return f"{self.chrom}_{self.pos}"
    
    def __post_init__(self) -> None:
        """Validate variant data after initialization."""
        if not self.chrom:
            raise ValueError("Chromosome cannot be empty")
        if self.pos <= 0:
            raise ValueError("Position must be positive")
        if not self.ref:
            raise ValueError("Reference allele cannot be empty")
        # Allow empty alt list for special cases (e.g., structural variants)
        if self.alt is None:
            self.alt = []


@dataclass
class ProcessingConfig:
    """Configuration for VCF processing.
    
    Attributes:
        vcf_file: Path to input VCF file
        output_file: Base path for output files
        target_id_file: Path to target IDs file (optional, if None processes all variants)
        miss_fmt: Format for missing genotypes (default: "NN")
        gt_sep: Separator for genotype alleles (default: "")
        threads: Number of processing threads (default: 4)
        batch_size: Batch size for processing (default: 10000)
        compress_output: Whether to compress output files (default: True)
        verbose: Enable verbose logging (default: False)
        quiet: Enable quiet mode (default: False)
        dry_run: Preview operations without executing (default: False)
        log_file: Optional log file path
        include_variant_type: Whether to include variant type annotation (default: False)
    """
    vcf_file: Path
    output_file: Path
    target_id_file: Optional[Path] = None
    miss_fmt: str = "./."
    gt_sep: str = ""
    threads: int = 4
    batch_size: int = 10000
    compress_output: bool = True
    verbose: bool = False
    quiet: bool = False
    dry_run: bool = False
    log_file: Optional[Path] = None
    include_variant_type: bool = False
    
    def __post_init__(self) -> None:
        """Validate configuration after initialization."""
        self.validate()
    
    def validate(self) -> None:
        """Validate configuration parameters.
        
        Raises:
            ValueError: If any parameter is invalid
            FileNotFoundError: If required input files don't exist
        """
        # Check input files exist
        if not self.vcf_file.exists():
            raise FileNotFoundError(f"VCF file not found: {self.vcf_file}")
        
        if self.target_id_file and not self.target_id_file.exists():
            raise FileNotFoundError(f"Target ID file not found: {self.target_id_file}")
        
        # Validate numeric parameters
        if self.threads <= 0:
            raise ValueError("Threads must be positive")
        
        if self.batch_size <= 0:
            raise ValueError("Batch size must be positive")
        
        # Validate miss_fmt
        if not self.miss_fmt:
            raise ValueError("Missing format cannot be empty")
        
        # Validate conflicting options
        if self.verbose and self.quiet:
            raise ValueError("Cannot enable both verbose and quiet modes")
        
        # Ensure output directory exists
        self.output_file.parent.mkdir(parents=True, exist_ok=True)
        
        logger.debug(f"Configuration validated: {self}")
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return {
            "vcf_file": str(self.vcf_file),
            "target_id_file": str(self.target_id_file) if self.target_id_file else None,
            "output_file": str(self.output_file),
            "miss_fmt": self.miss_fmt,
            "gt_sep": self.gt_sep,
            "threads": self.threads,
            "batch_size": self.batch_size,
            "compress_output": self.compress_output,
            "verbose": self.verbose,
            "quiet": self.quiet,
            "dry_run": self.dry_run,
            "log_file": str(self.log_file) if self.log_file else None,
            "include_variant_type": self.include_variant_type,
        }


@dataclass
class ProcessingResult:
    """Result of VCF processing operation.
    
    Attributes:
        total_variants: Total number of variants in input
        processed_variants: Number of variants successfully processed
        skipped_variants: Number of variants skipped
        errors: List of error messages encountered
        processing_time: Total processing time in seconds
        output_files: List of generated output files
        config: Configuration used for processing
    """
    total_variants: int = 0
    processed_variants: int = 0
    skipped_variants: int = 0
    errors: List[str] = field(default_factory=list)
    processing_time: float = 0.0
    output_files: List[Path] = field(default_factory=list)
    config: Optional[ProcessingConfig] = None
    
    @property
    def success_rate(self) -> float:
        """Calculate processing success rate."""
        if self.total_variants == 0:
            return 0.0
        return self.processed_variants / self.total_variants
    
    @property
    def has_errors(self) -> bool:
        """Check if any errors occurred."""
        return len(self.errors) > 0
    
    def add_error(self, error: str) -> None:
        """Add an error message."""
        self.errors.append(error)
        logger.error(error)
    
    def summary(self) -> str:
        """Generate a summary report of processing results."""
        lines = [
            "=== VCF Processing Summary ===",
            f"Total variants: {self.total_variants:,}",
            f"Processed: {self.processed_variants:,}",
            f"Skipped: {self.skipped_variants:,}",
            f"Success rate: {self.success_rate:.1%}",
            f"Processing time: {self.processing_time:.2f}s",
        ]
        
        if self.output_files:
            lines.append("Output files:")
            for file_path in self.output_files:
                size = file_path.stat().st_size if file_path.exists() else 0
                lines.append(f"  - {file_path} ({size:,} bytes)")
        
        if self.has_errors:
            lines.append(f"Errors ({len(self.errors)}):")
            for error in self.errors[:5]:  # Show first 5 errors
                lines.append(f"  - {error}")
            if len(self.errors) > 5:
                lines.append(f"  ... and {len(self.errors) - 5} more errors")
        
        return "\n".join(lines)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert result to dictionary."""
        return {
            "total_variants": self.total_variants,
            "processed_variants": self.processed_variants,
            "skipped_variants": self.skipped_variants,
            "success_rate": self.success_rate,
            "processing_time": self.processing_time,
            "errors": self.errors,
            "output_files": [str(f) for f in self.output_files],
            "has_errors": self.has_errors,
        }


class ProcessingTimer:
    """Context manager for timing processing operations."""
    
    def __init__(self, result: ProcessingResult) -> None:
        """Initialize timer with result object to update."""
        self.result = result
        self.start_time: float = 0.0
    
    def __enter__(self) -> "ProcessingTimer":
        """Start timing."""
        self.start_time = time.time()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Stop timing and update result."""
        self.result.processing_time = time.time() - self.start_time