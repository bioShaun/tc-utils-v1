"""Variant filtering functionality."""

from pathlib import Path
from typing import Set, Optional, List

from loguru import logger

from .config import VariantInfo


class VariantFilter:
    """Handles filtering of variants based on target IDs and other criteria.
    
    This class provides methods to load target variant IDs from files
    and filter variants based on various criteria.
    """
    
    def __init__(self, target_ids: Optional[Set[str]] = None):
        """Initialize variant filter.
        
        Args:
            target_ids: Optional set of target variant IDs to filter by
        """
        self.target_ids = target_ids or set()
        self._stats = {
            'total_checked': 0,
            'passed_filter': 0,
            'failed_filter': 0
        }
    
    @classmethod
    def from_file(cls, target_file: Optional[Path]) -> 'VariantFilter':
        """Create VariantFilter from target ID file.
        
        Args:
            target_file: Path to file containing variant IDs (one per line), or None to process all variants
            
        Returns:
            VariantFilter instance with loaded target IDs
            
        Raises:
            FileNotFoundError: If target file doesn't exist
            ValueError: If target file is empty or malformed
        """
        target_ids = load_target_ids(target_file)
        return cls(target_ids)
    
    def add_target_id(self, variant_id: str) -> None:
        """Add a single target ID to the filter.
        
        Args:
            variant_id: Variant ID to add
        """
        self.target_ids.add(variant_id)
    
    def add_target_ids(self, variant_ids: Set[str]) -> None:
        """Add multiple target IDs to the filter.
        
        Args:
            variant_ids: Set of variant IDs to add
        """
        self.target_ids.update(variant_ids)
    
    def remove_target_id(self, variant_id: str) -> None:
        """Remove a target ID from the filter.
        
        Args:
            variant_id: Variant ID to remove
        """
        self.target_ids.discard(variant_id)
    
    def clear_target_ids(self) -> None:
        """Clear all target IDs from the filter."""
        self.target_ids.clear()
    
    def passes_filter(self, variant: VariantInfo) -> bool:
        """Check if a variant passes the filter criteria.
        
        Args:
            variant: VariantInfo object to check
            
        Returns:
            True if variant passes filter, False otherwise
        """
        self._stats['total_checked'] += 1
        
        # If no target IDs specified, all variants pass
        if not self.target_ids:
            self._stats['passed_filter'] += 1
            return True
        
        # Check if variant ID is in target set
        if variant.variant_id in self.target_ids:
            self._stats['passed_filter'] += 1
            return True
        
        self._stats['failed_filter'] += 1
        return False
    
    def filter_variants(self, variants: List[VariantInfo]) -> List[VariantInfo]:
        """Filter a list of variants.
        
        Args:
            variants: List of VariantInfo objects to filter
            
        Returns:
            List of variants that pass the filter
        """
        return [variant for variant in variants if self.passes_filter(variant)]
    
    def get_target_count(self) -> int:
        """Get number of target IDs in the filter.
        
        Returns:
            Number of target IDs
        """
        return len(self.target_ids)
    
    def get_stats(self) -> dict:
        """Get filtering statistics.
        
        Returns:
            Dictionary with filtering statistics
        """
        return self._stats.copy()
    
    def reset_stats(self) -> None:
        """Reset filtering statistics."""
        self._stats = {
            'total_checked': 0,
            'passed_filter': 0,
            'failed_filter': 0
        }
    
    def validate_target_ids(self) -> List[str]:
        """Validate target IDs format.
        
        Returns:
            List of invalid target IDs (empty if all valid)
        """
        invalid_ids = []
        
        for target_id in self.target_ids:
            # Check basic format: should contain at least one underscore
            if '_' not in target_id:
                invalid_ids.append(target_id)
                continue
            
            # Check if it follows CHROM_POS pattern
            parts = target_id.split('_')
            if len(parts) != 2:
                invalid_ids.append(target_id)
                continue
            
            chrom, pos_str = parts
            if not chrom or not pos_str:
                invalid_ids.append(target_id)
                continue
            
            # Check if position is numeric
            try:
                pos = int(pos_str)
                if pos <= 0:
                    invalid_ids.append(target_id)
            except ValueError:
                invalid_ids.append(target_id)
        
        if invalid_ids:
            logger.warning(f"Found {len(invalid_ids)} invalid target IDs: {invalid_ids[:5]}...")
        
        return invalid_ids


def load_target_ids(target_file: Optional[Path]) -> Set[str]:
    """Load target variant IDs from file.
    
    Args:
        target_file: Path to file containing variant IDs (one per line), or None to process all variants
        
    Returns:
        Set of variant IDs, empty set if target_file is None
        
    Raises:
        FileNotFoundError: If target file doesn't exist
        ValueError: If target file is empty or malformed
    """
    if target_file is None:
        logger.info("No target file specified, will process all variants")
        return set()
        
    if not target_file.exists():
        raise FileNotFoundError(f"Target ID file not found: {target_file}")
    
    target_ids = set()
    invalid_lines = []
    
    with open(target_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            
            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue
            
            # Basic validation: should not contain whitespace
            if ' ' in line or '\t' in line:
                invalid_lines.append((line_num, line))
                continue
            
            target_ids.add(line)
    
    if invalid_lines:
        logger.warning(f"Skipped {len(invalid_lines)} invalid lines in {target_file}")
        for line_num, line in invalid_lines[:5]:  # Show first 5 invalid lines
            logger.warning(f"  Line {line_num}: {repr(line)}")
    
    if not target_ids:
        raise ValueError(f"No valid target IDs found in {target_file}")
    
    logger.info(f"Loaded {len(target_ids)} target IDs from {target_file}")
    return target_ids


def validate_target_file(target_file: Path) -> bool:
    """Validate target ID file format.
    
    Args:
        target_file: Path to target ID file
        
    Returns:
        True if file is valid, False otherwise
    """
    try:
        target_ids = load_target_ids(target_file)
        
        # Additional validation using VariantFilter
        filter_obj = VariantFilter(target_ids)
        invalid_ids = filter_obj.validate_target_ids()
        
        if invalid_ids:
            logger.error(f"Target file contains {len(invalid_ids)} invalid IDs")
            return False
        
        return True
        
    except Exception as e:
        logger.error(f"Target file validation failed: {e}")
        return False