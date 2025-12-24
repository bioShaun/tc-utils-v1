"""Variant transformation functionality."""

from typing import List, Dict, Any
from loguru import logger

from .config import VariantInfo


class VariantTransformer:
    """Handles transformation of VCF variants to annotation format.
    
    This class provides methods to transform VCF variant representations
    into more readable annotation formats, handling SNPs, insertions,
    deletions, and complex variants.
    """
    
    def __init__(self, strict_mode: bool = True):
        """Initialize variant transformer.
        
        Args:
            strict_mode: If True, raises exceptions for complex variants.
                        If False, returns original representation for complex variants.
        """
        self.strict_mode = strict_mode
        self._stats = {
            'snps': 0,
            'insertions': 0,
            'deletions': 0,
            'mnps': 0,
            'complex': 0,
            'errors': 0
        }
    
    def transform_one_alt(self, ref: str, alt: str) -> str:
        """Transform a single VCF variant to annotation format.
        
        Args:
            ref: Reference allele
            alt: Alternative allele
            
        Returns:
            Transformed variant representation
            
        Raises:
            ValueError: If ref/alt are empty or complex variant in strict mode
        """
        # Handle empty inputs
        if not ref or not alt:
            if self.strict_mode:
                raise ValueError("ref and alt cannot be empty")
            else:
                self._stats['errors'] += 1
                logger.error(f"Empty input: ref='{ref}', alt='{alt}'")
                return alt
        
        # Handle symbolic alleles
        if alt == "*":
            self._stats['deletions'] += 1
            return f"del{ref}"
        
        try:
            # Remove common prefix
            prefix_len = 0
            min_len = min(len(ref), len(alt))
            while prefix_len < min_len and ref[prefix_len] == alt[prefix_len]:
                prefix_len += 1
            
            # Remove common suffix (but don't overlap with prefix removal)
            suffix_len = 0
            # Only look at the remaining parts after prefix removal
            ref_remaining = ref[prefix_len:]
            alt_remaining = alt[prefix_len:]
            min_remaining = min(len(ref_remaining), len(alt_remaining))
            
            while (suffix_len < min_remaining and 
                   ref_remaining[-(suffix_len + 1)] == alt_remaining[-(suffix_len + 1)]):
                suffix_len += 1
            
            # Extract changed parts
            if suffix_len > 0:
                ref_changed = ref_remaining[:-suffix_len]
                alt_changed = alt_remaining[:-suffix_len]
            else:
                ref_changed = ref_remaining
                alt_changed = alt_remaining
            
            # Pure deletion
            if ref_changed and not alt_changed:
                self._stats['deletions'] += 1
                return f"del{ref_changed}"
            
            # Pure insertion
            if not ref_changed and alt_changed:
                self._stats['insertions'] += 1
                return f"ins{alt_changed}"
            
            # SNP (single nucleotide substitution)
            if len(ref_changed) == 1 and len(alt_changed) == 1:
                self._stats['snps'] += 1
                return alt_changed
            
            # MNP (multi-nucleotide substitution, same length)
            if len(ref_changed) == len(alt_changed) and ref_changed and alt_changed:
                self._stats['mnps'] += 1
                return alt_changed
            
            # Complex variant (delins) - different lengths and both non-empty
            if ref_changed and alt_changed and len(ref_changed) != len(alt_changed):
                self._stats['complex'] += 1
                if self.strict_mode:
                    raise ValueError(
                        f"Complex variant (delins) detected: "
                        f"ref={ref}, alt={alt}, "
                        f"ref_changed={ref_changed}, alt_changed={alt_changed}"
                    )
                else:
                    logger.warning(f"Complex variant detected, returning original: {ref}->{alt}")
                    return alt
            
            # Fallback - return alt if we can't classify
            self._stats['complex'] += 1
            if self.strict_mode:
                raise ValueError(
                    f"Unclassified variant: ref={ref}, alt={alt}, "
                    f"ref_changed={ref_changed}, alt_changed={alt_changed}"
                )
            else:
                logger.warning(f"Unclassified variant, returning original: {ref}->{alt}")
                return alt
                
        except Exception as e:
            self._stats['errors'] += 1
            if self.strict_mode:
                raise
            else:
                logger.error(f"Error transforming variant {ref}->{alt}: {e}")
                return alt
    
    def transform_alt_alleles(self, ref: str, alt_alleles: List[str]) -> List[str]:
        """Transform multiple alternative alleles.
        
        Args:
            ref: Reference allele
            alt_alleles: List of alternative alleles (can be empty)
            
        Returns:
            List of transformed alternative alleles
        """
        if not alt_alleles:
            return []
        return [self.transform_one_alt(ref, alt) for alt in alt_alleles]
    
    def transform_alt_string(self, ref: str, alt_string: str) -> str:
        """Transform comma-separated alternative alleles string.
        
        Args:
            ref: Reference allele
            alt_string: Comma-separated alternative alleles
            
        Returns:
            Comma-separated transformed alternative alleles
        """
        alt_list = alt_string.split(",")
        transformed_list = self.transform_alt_alleles(ref, alt_list)
        return ",".join(transformed_list)
    
    def transform_variant(self, variant: VariantInfo) -> VariantInfo:
        """Transform a VariantInfo object.
        
        Args:
            variant: VariantInfo object to transform
            
        Returns:
            New VariantInfo object with transformed alleles
        """
        # Transform alternative alleles
        transformed_alt = self.transform_alt_alleles(variant.ref, variant.alt)
        
        # Transform reference allele (take first base for multi-base refs)
        transformed_ref = variant.ref[0] if len(variant.ref) > 1 else variant.ref
        
        # Create new variant with transformed alleles
        return VariantInfo(
            chrom=variant.chrom,
            pos=variant.pos,
            ref=transformed_ref,
            alt=transformed_alt,
            variant_id=variant.variant_id,
            genotypes=variant.genotypes
        )
    
    def get_stats(self) -> Dict[str, int]:
        """Get transformation statistics.
        
        Returns:
            Dictionary with transformation statistics
        """
        return self._stats.copy()
    
    def reset_stats(self) -> None:
        """Reset transformation statistics."""
        self._stats = {
            'snps': 0,
            'insertions': 0,
            'deletions': 0,
            'mnps': 0,
            'complex': 0,
            'errors': 0
        }
    
    def get_summary(self) -> str:
        """Get a summary of transformation statistics.
        
        Returns:
            Human-readable summary string
        """
        stats = self.get_stats()
        total = sum(stats.values())
        
        if total == 0:
            return "No variants processed"
        
        lines = [
            f"Transformation Summary (Total: {total})",
            f"  SNPs: {stats['snps']} ({stats['snps']/total*100:.1f}%)",
            f"  Insertions: {stats['insertions']} ({stats['insertions']/total*100:.1f}%)",
            f"  Deletions: {stats['deletions']} ({stats['deletions']/total*100:.1f}%)",
            f"  MNPs: {stats['mnps']} ({stats['mnps']/total*100:.1f}%)",
            f"  Complex: {stats['complex']} ({stats['complex']/total*100:.1f}%)",
            f"  Errors: {stats['errors']} ({stats['errors']/total*100:.1f}%)"
        ]
        
        return "\n".join(lines)
    
    @staticmethod
    def classify_variant_type(ref: str, alt: str) -> str:
        """Classify a single variant type based on REF and ALT.
        
        Args:
            ref: Reference allele
            alt: Alternative allele
            
        Returns:
            Variant type: "SNP", "INDEL", "MNP", "REF", or "UNKNOWN"
        """
        # Handle empty inputs
        if not ref or not alt:
            return "UNKNOWN"
        
        # Handle symbolic alleles
        if alt == "*":
            return "INDEL"  # Deletion
        
        # Normalize case for comparison
        ref_upper = ref.upper()
        alt_upper = alt.upper()
        
        # Reference allele (no variation)
        if ref_upper == alt_upper:
            return "REF"
        
        # Single nucleotide polymorphism
        if len(ref) == 1 and len(alt) == 1:
            return "SNP"
        
        # Insertion or deletion (different lengths)
        if len(ref) != len(alt):
            return "INDEL"
        
        # Multi-nucleotide polymorphism (same length, >1 base)
        if len(ref) == len(alt) and len(ref) > 1:
            return "MNP"
        
        # Fallback for edge cases
        return "UNKNOWN"
    
    @staticmethod
    def classify_multi_allelic_types(ref: str, alts: List[str]) -> str:
        """Classify multi-allelic variant types.
        
        Args:
            ref: Reference allele
            alts: List of alternative alleles
            
        Returns:
            Combined variant types separated by pipe (|) and sorted as SNP|INDEL|MNP|REF
        """
        if not alts:
            return "REF"
        
        types = set()
        for alt in alts:
            variant_type = VariantTransformer.classify_variant_type(ref, alt)
            if variant_type != "UNKNOWN":  # Only include known types
                types.add(variant_type)
        
        # Sort types in predefined order
        type_order = ["SNP", "INDEL", "MNP", "REF"]
        sorted_types = [t for t in type_order if t in types]
        
        # If no known types found, return UNKNOWN
        if not sorted_types:
            return "UNKNOWN"
        
        return "|".join(sorted_types)
    
    def get_variant_type_column(self, variants: List[VariantInfo]) -> List[str]:
        """Get variant type column for a list of variants.
        
        Args:
            variants: List of VariantInfo objects
            
        Returns:
            List of variant type strings for each variant
        """
        variant_types = []
        for variant in variants:
            variant_type = self.classify_multi_allelic_types(variant.ref, variant.alt)
            variant_types.append(variant_type)
        return variant_types


def validate_variant_alleles(ref: str, alt_alleles: List[str]) -> List[str]:
    """Validate variant alleles and return list of issues.
    
    Args:
        ref: Reference allele
        alt_alleles: List of alternative alleles
        
    Returns:
        List of validation issues (empty if all valid)
    """
    issues = []
    
    # Check reference allele
    if not ref:
        issues.append("Reference allele is empty")
    elif not all(c in 'ATCGN' for c in ref.upper()):
        issues.append(f"Reference allele contains invalid characters: {ref}")
    
    # Check alternative alleles
    if not alt_alleles:
        issues.append("No alternative alleles provided")
    else:
        for i, alt in enumerate(alt_alleles):
            if not alt:
                issues.append(f"Alternative allele {i+1} is empty")
            elif alt != "*" and not all(c in 'ATCGN' for c in alt.upper()):
                issues.append(f"Alternative allele {i+1} contains invalid characters: {alt}")
    
    return issues