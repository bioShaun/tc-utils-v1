"""VCF file reading with cyvcf2."""

from pathlib import Path
from typing import Iterator, List, Optional, Set

from cyvcf2 import VCF
from loguru import logger

from .config import VariantInfo


class VCFReader:
    """VCF file reader using cyvcf2 for efficient VCF parsing.
    
    This class provides methods to read VCF files, extract sample names,
    and iterate over variants with optional filtering by target IDs.
    """
    
    def __init__(self, vcf_path: Path, threads: int = 1):
        """Initialize VCF reader.
        
        Args:
            vcf_path: Path to VCF file (can be compressed)
            threads: Number of threads for reading (cyvcf2 parameter)
        """
        self.vcf_path = vcf_path
        self.threads = threads
        self._vcf = None
        self._sample_names = None
        
        if not vcf_path.exists():
            raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    
    def __enter__(self):
        """Context manager entry."""
        self._vcf = VCF(str(self.vcf_path), threads=self.threads)
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        if self._vcf:
            self._vcf.close()
            self._vcf = None
    
    @property
    def sample_names(self) -> List[str]:
        """Get sample names from VCF header.
        
        Returns:
            List of sample names in the order they appear in the VCF
        """
        if self._sample_names is None:
            if self._vcf is None:
                # Temporarily open VCF to get sample names
                vcf = VCF(str(self.vcf_path))
                self._sample_names = list(vcf.samples)
                vcf.close()
            else:
                self._sample_names = list(self._vcf.samples)
        return self._sample_names
    
    def get_sample_count(self) -> int:
        """Get number of samples in VCF.
        
        Returns:
            Number of samples
        """
        return len(self.sample_names)
    
    def iter_variants(self, target_ids: Optional[Set[str]] = None) -> Iterator[VariantInfo]:
        """Iterate over variants in the VCF file.
        
        Args:
            target_ids: Optional set of variant IDs to filter by.
                       If None, all variants are yielded.
        
        Yields:
            VariantInfo objects for each variant
        """
        if self._vcf is None:
            raise RuntimeError("VCFReader must be used as context manager")
        
        for variant in self._vcf:
            # Generate variant ID in format CHROM_POS
            variant_id = f"{variant.CHROM}_{variant.POS}"
            
            # Filter by target IDs if provided
            if target_ids is not None and variant_id not in target_ids:
                continue
            
            # Extract genotypes for all samples
            genotypes = []
            for gt in variant.genotypes:
                # Convert cyvcf2 genotype format to VCF format
                if gt[0] == -1 or gt[1] == -1:
                    # Missing genotype
                    genotypes.append("./.")
                else:
                    # Format as allele1/allele2, normalize phasing to unphased
                    gt_str = f"{gt[0]}/{gt[1]}"
                    genotypes.append(gt_str)
            
            # Handle multi-allelic variants
            alt_alleles = variant.ALT if isinstance(variant.ALT, list) else [variant.ALT]
            
            yield VariantInfo(
                chrom=variant.CHROM,
                pos=variant.POS,
                ref=variant.REF,
                alt=alt_alleles,
                variant_id=variant_id,
                genotypes=genotypes
            )
    
    def count_variants(self, target_ids: Optional[Set[str]] = None) -> int:
        """Count variants in the VCF file.
        
        Args:
            target_ids: Optional set of variant IDs to filter by
            
        Returns:
            Number of variants (filtered if target_ids provided)
        """
        count = 0
        for _ in self.iter_variants(target_ids):
            count += 1
        return count
    
    def validate_vcf(self) -> bool:
        """Validate that the VCF file can be opened and read.
        
        Returns:
            True if VCF is valid, False otherwise
        """
        try:
            with VCF(str(self.vcf_path)) as vcf:
                # Try to read first variant
                for _ in vcf:
                    break
            return True
        except Exception as e:
            logger.error(f"VCF validation failed: {e}")
            return False


def load_target_ids(target_file: Path) -> Set[str]:
    """Load target variant IDs from file.
    
    Args:
        target_file: Path to file containing variant IDs (one per line)
        
    Returns:
        Set of variant IDs
        
    Raises:
        FileNotFoundError: If target file doesn't exist
        ValueError: If target file is empty or malformed
    """
    if not target_file.exists():
        raise FileNotFoundError(f"Target ID file not found: {target_file}")
    
    target_ids = set()
    with open(target_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if line and not line.startswith('#'):
                target_ids.add(line)
    
    if not target_ids:
        raise ValueError(f"No valid target IDs found in {target_file}")
    
    logger.info(f"Loaded {len(target_ids)} target IDs from {target_file}")
    return target_ids