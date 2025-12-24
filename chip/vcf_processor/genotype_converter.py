"""Genotype conversion functionality."""

from typing import List, Dict, Any, Optional, Iterator
import pandas as pd
from loguru import logger

from .config import VariantInfo


class GenotypeConverter:
    """Handles conversion of VCF genotypes to sequence format.
    
    This class provides methods to convert VCF genotype calls (like 0/1, 1/1)
    to actual sequence representations using reference and alternative alleles.
    """
    
    def __init__(self, miss_fmt: str = "NN", gt_sep: str = ""):
        """Initialize genotype converter.
        
        Args:
            miss_fmt: Format for missing genotypes (default: "NN")
            gt_sep: Separator for heterozygous genotypes (default: "")
        """
        self.miss_fmt = miss_fmt
        self.gt_sep = gt_sep
        self._stats = {
            'total_genotypes': 0,
            'missing_genotypes': 0,
            'homozygous_ref': 0,
            'homozygous_alt': 0,
            'heterozygous': 0,
            'errors': 0
        }
    
    def convert_genotype(self, genotype: str, ref: str, alt_alleles: List[str]) -> str:
        """Convert a single genotype to sequence format.
        
        Args:
            genotype: VCF genotype string (e.g., "0/1", "1/1", "./.")
            ref: Reference allele
            alt_alleles: List of alternative alleles
            
        Returns:
            Converted genotype sequence
        """
        self._stats['total_genotypes'] += 1
        
        try:
            # Handle missing genotype
            if genotype == "./." or genotype == ".|.":
                self._stats['missing_genotypes'] += 1
                return self.miss_fmt
            
            # Parse genotype
            if "/" in genotype:
                allele1_str, allele2_str = genotype.split("/")
            elif "|" in genotype:
                allele1_str, allele2_str = genotype.split("|")
            else:
                # Assume homozygous if no separator
                allele1_str = allele2_str = genotype
            
            # Build allele list (0=ref, 1=alt1, 2=alt2, etc.)
            allele_list = [ref] + alt_alleles
            
            # Convert allele indices to sequences
            allele1_seq = "N" if allele1_str == "." else allele_list[int(allele1_str)]
            allele2_seq = "N" if allele2_str == "." else allele_list[int(allele2_str)]
            
            # Format output based on genotype type
            if allele1_str == allele2_str:
                # Homozygous
                if allele1_str == "0":
                    self._stats['homozygous_ref'] += 1
                else:
                    self._stats['homozygous_alt'] += 1
                
                # For multi-base alleles, return just the sequence
                if len(allele1_seq) > 1:
                    return allele1_seq
                return f"{allele1_seq}{self.gt_sep}{allele2_seq}"
            else:
                # Heterozygous
                self._stats['heterozygous'] += 1
                
                # For multi-base alleles, use "/" separator
                if len(allele1_seq) > 1 or len(allele2_seq) > 1:
                    return f"{allele1_seq}/{allele2_seq}"
                return f"{allele1_seq}{self.gt_sep}{allele2_seq}"
                
        except (ValueError, IndexError) as e:
            self._stats['errors'] += 1
            logger.error(f"Error converting genotype {genotype}: {e}")
            return self.miss_fmt
    
    def convert_variant_genotypes(self, variant: VariantInfo, sample_names: List[str]) -> Dict[str, str]:
        """Convert all genotypes for a variant.
        
        Args:
            variant: VariantInfo object containing genotype data
            sample_names: List of sample names corresponding to genotypes
            
        Returns:
            Dictionary mapping sample names to converted genotypes
        """
        if len(variant.genotypes) != len(sample_names):
            raise ValueError(f"Genotype count ({len(variant.genotypes)}) doesn't match sample count ({len(sample_names)})")
        
        result = {}
        for sample_name, genotype in zip(sample_names, variant.genotypes):
            converted = self.convert_genotype(genotype, variant.ref, variant.alt)
            result[sample_name] = converted
        
        return result
    
    def convert_batch(self, variants: List[VariantInfo], sample_names: List[str]) -> pd.DataFrame:
        """Convert a batch of variants to DataFrame format.
        
        Args:
            variants: List of VariantInfo objects
            sample_names: List of sample names
            
        Returns:
            DataFrame with converted genotypes
        """
        rows = []
        
        for variant in variants:
            # Convert genotypes for this variant
            converted_gts = self.convert_variant_genotypes(variant, sample_names)
            
            # Create row with location info and genotypes
            row = {
                'CHROM': variant.chrom,
                'POS': variant.pos,
                'REF': variant.ref,
                'ALT': ','.join(variant.alt) if isinstance(variant.alt, list) else variant.alt
            }
            row.update(converted_gts)
            rows.append(row)
        
        return pd.DataFrame(rows)
    
    def process_dataframe(self, df: pd.DataFrame, location_cols: List[str] = None) -> pd.DataFrame:
        """Process a DataFrame with VCF-style genotype data.
        
        Args:
            df: DataFrame with CHROM, POS, REF, ALT columns and sample genotype columns
            location_cols: List of location column names (default: ['CHROM', 'POS', 'REF', 'ALT'])
            
        Returns:
            DataFrame with converted genotypes
        """
        if location_cols is None:
            location_cols = ['CHROM', 'POS', 'REF', 'ALT']
        
        # Make a copy to avoid modifying original
        result_df = df.copy()
        
        # Get sample columns (all columns except location columns)
        sample_cols = [col for col in df.columns if col not in location_cols]
        
        # Convert each genotype
        for _, row in df.iterrows():
            ref = row['REF']
            alt_alleles = row['ALT'].split(',') if isinstance(row['ALT'], str) else [row['ALT']]
            
            for sample_col in sample_cols:
                genotype = row[sample_col]
                converted = self.convert_genotype(genotype, ref, alt_alleles)
                result_df.at[row.name, sample_col] = converted
        
        return result_df
    
    def get_stats(self) -> Dict[str, int]:
        """Get conversion statistics.
        
        Returns:
            Dictionary with conversion statistics
        """
        return self._stats.copy()
    
    def reset_stats(self) -> None:
        """Reset conversion statistics."""
        self._stats = {
            'total_genotypes': 0,
            'missing_genotypes': 0,
            'homozygous_ref': 0,
            'homozygous_alt': 0,
            'heterozygous': 0,
            'errors': 0
        }
    
    def get_summary(self) -> str:
        """Get a summary of conversion statistics.
        
        Returns:
            Human-readable summary string
        """
        stats = self.get_stats()
        total = stats['total_genotypes']
        
        if total == 0:
            return "No genotypes processed"
        
        lines = [
            f"Genotype Conversion Summary (Total: {total})",
            f"  Missing: {stats['missing_genotypes']} ({stats['missing_genotypes']/total*100:.1f}%)",
            f"  Homozygous Ref: {stats['homozygous_ref']} ({stats['homozygous_ref']/total*100:.1f}%)",
            f"  Homozygous Alt: {stats['homozygous_alt']} ({stats['homozygous_alt']/total*100:.1f}%)",
            f"  Heterozygous: {stats['heterozygous']} ({stats['heterozygous']/total*100:.1f}%)",
            f"  Errors: {stats['errors']} ({stats['errors']/total*100:.1f}%)"
        ]
        
        return "\n".join(lines)


class BatchProcessor:
    """Handles batch processing of large VCF datasets with memory management."""
    
    def __init__(self, batch_size: int = 1000, progress_callback: Optional[callable] = None):
        """Initialize batch processor.
        
        Args:
            batch_size: Number of variants to process in each batch
            progress_callback: Optional callback function for progress reporting
        """
        self.batch_size = batch_size
        self.progress_callback = progress_callback
        self._processed_count = 0
        self._total_count = 0
    
    def process_variants(self, variants: Iterator[VariantInfo], 
                        sample_names: List[str],
                        converter: GenotypeConverter) -> Iterator[pd.DataFrame]:
        """Process variants in batches.
        
        Args:
            variants: Iterator of VariantInfo objects
            sample_names: List of sample names
            converter: GenotypeConverter instance
            
        Yields:
            DataFrames with converted genotypes for each batch
        """
        batch = []
        
        for variant in variants:
            batch.append(variant)
            
            if len(batch) >= self.batch_size:
                # Process batch
                df = converter.convert_batch(batch, sample_names)
                self._processed_count += len(batch)
                
                if self.progress_callback:
                    self.progress_callback(self._processed_count, len(batch))
                
                yield df
                batch = []
        
        # Process remaining variants
        if batch:
            df = converter.convert_batch(batch, sample_names)
            self._processed_count += len(batch)
            
            if self.progress_callback:
                self.progress_callback(self._processed_count, len(batch))
            
            yield df
    
    def get_processed_count(self) -> int:
        """Get number of variants processed so far.
        
        Returns:
            Number of processed variants
        """
        return self._processed_count
    
    def reset_count(self) -> None:
        """Reset processed count."""
        self._processed_count = 0