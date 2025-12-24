"""Output writing functionality for VCF processor."""

import gzip
from pathlib import Path
from typing import Any, Dict, List, Optional, TextIO, Union

import pandas as pd
from loguru import logger

from .config import ProcessingConfig, ProcessingResult


class OutputWriter:
    """Handles output file writing with compression and batch support.
    
    This class manages writing genotype and sequence data to output files,
    supporting both compressed and uncompressed formats, proper file naming
    conventions, and append mode for batch processing.
    
    Attributes:
        config: Processing configuration
        result: Processing result to track output files
        _gt_file: File handle for genotype output
        _seq_file: File handle for sequence output
        _is_first_batch: Flag to track if this is the first batch
    """
    
    def __init__(self, config: ProcessingConfig, result: ProcessingResult) -> None:
        """Initialize OutputWriter.
        
        Args:
            config: Processing configuration
            result: Processing result to track output files
        """
        self.config = config
        self.result = result
        self._gt_file: Optional[TextIO] = None
        self._seq_file: Optional[TextIO] = None
        self._is_first_batch = True
        
        # Generate output file paths with professional suffixes
        self._gt_path = self._get_output_path("genotype_codes")  # .genotype_codes.tsv for VCF format genotypes (0/0, 0/1, etc.)
        self._seq_path = self._get_output_path("genotype_bases")  # .genotype_bases.tsv for sequence format (ATCG)
        
        logger.debug(f"OutputWriter initialized with paths: {self._gt_path}, {self._seq_path}")
    
    def _get_output_path(self, suffix: str) -> Path:
        """Generate output file path with proper naming convention.
        
        Args:
            suffix: File suffix (e.g., 'genotype_codes', 'genotype_bases')
            
        Returns:
            Path to output file with proper extension
        """
        base_path = self.config.output_file
        extension = ".tsv.gz" if self.config.compress_output else ".tsv"
        return Path(f"{base_path}.{suffix}{extension}")
    
    def __enter__(self) -> "OutputWriter":
        """Enter context manager and open output files."""
        try:
            # Determine file mode
            mode = "wt" if self._is_first_batch else "at"
            
            # Open files with appropriate compression
            if self.config.compress_output:
                self._gt_file = gzip.open(self._gt_path, mode, encoding="utf-8")
                self._seq_file = gzip.open(self._seq_path, mode, encoding="utf-8")
            else:
                self._gt_file = open(self._gt_path, mode, encoding="utf-8")
                self._seq_file = open(self._seq_path, mode, encoding="utf-8")
            
            logger.debug(f"Opened output files in mode '{mode}'")
            return self
            
        except Exception as e:
            error_msg = f"Failed to open output files: {e}"
            self.result.add_error(error_msg)
            raise
    
    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Exit context manager and close output files."""
        try:
            if self._gt_file:
                self._gt_file.close()
                self._gt_file = None
            
            if self._seq_file:
                self._seq_file.close()
                self._seq_file = None
            
            # Add output files to result if they exist and have content
            if self._gt_path.exists() and self._gt_path.stat().st_size > 0:
                if self._gt_path not in self.result.output_files:
                    self.result.output_files.append(self._gt_path)
            
            if self._seq_path.exists() and self._seq_path.stat().st_size > 0:
                if self._seq_path not in self.result.output_files:
                    self.result.output_files.append(self._seq_path)
            
            logger.debug("Closed output files")
            
        except Exception as e:
            error_msg = f"Error closing output files: {e}"
            self.result.add_error(error_msg)
            logger.error(error_msg)
    
    def write_batch(self, gt_df: pd.DataFrame, seq_df: pd.DataFrame) -> None:
        """Write a batch of genotype and sequence data.
        
        Args:
            gt_df: DataFrame containing genotype data
            seq_df: DataFrame containing sequence data
            
        Raises:
            RuntimeError: If files are not open or writing fails
        """
        if not self._gt_file or not self._seq_file:
            raise RuntimeError("Output files not open. Use within context manager.")
        
        try:
            # Write headers only for first batch
            include_header = self._is_first_batch
            
            # Write genotype data
            gt_df.to_csv(
                self._gt_file,
                sep="\t",
                index=False,
                header=include_header,
                lineterminator="\n"
            )
            
            # Write sequence data
            seq_df.to_csv(
                self._seq_file,
                sep="\t",
                index=False,
                header=include_header,
                lineterminator="\n"
            )
            
            # Update batch tracking
            self._is_first_batch = False
            
            logger.debug(f"Wrote batch: {len(gt_df)} genotype rows, {len(seq_df)} sequence rows")
            
        except Exception as e:
            error_msg = f"Failed to write batch: {e}"
            self.result.add_error(error_msg)
            raise RuntimeError(error_msg) from e
    
    def write_dataframes(self, gt_df: pd.DataFrame, seq_df: pd.DataFrame) -> None:
        """Write complete DataFrames to output files.
        
        This is a convenience method for writing complete datasets in one call.
        
        Args:
            gt_df: DataFrame containing genotype data
            seq_df: DataFrame containing sequence data
        """
        with self:
            self.write_batch(gt_df, seq_df)
    
    def get_output_paths(self) -> Dict[str, Path]:
        """Get dictionary of output file paths.
        
        Returns:
            Dictionary mapping file type to path
        """
        return {
            "genotype_codes": self._gt_path,
            "genotype_bases": self._seq_path,
        }
    
    def validate_output_files(self) -> bool:
        """Validate that output files exist and have content.
        
        Returns:
            True if all output files are valid, False otherwise
        """
        try:
            for file_path in [self._gt_path, self._seq_path]:
                if not file_path.exists():
                    self.result.add_error(f"Output file does not exist: {file_path}")
                    return False
                
                if file_path.stat().st_size == 0:
                    self.result.add_error(f"Output file is empty: {file_path}")
                    return False
            
            logger.debug("Output files validation passed")
            return True
            
        except Exception as e:
            error_msg = f"Error validating output files: {e}"
            self.result.add_error(error_msg)
            return False
    
    def cleanup_on_error(self) -> None:
        """Clean up output files if processing fails.
        
        This method removes incomplete output files to prevent
        confusion about processing status.
        """
        try:
            for file_path in [self._gt_path, self._seq_path]:
                if file_path.exists():
                    file_path.unlink()
                    logger.debug(f"Cleaned up incomplete file: {file_path}")
        
        except Exception as e:
            logger.warning(f"Failed to cleanup files: {e}")


class BatchOutputWriter(OutputWriter):
    """OutputWriter optimized for batch processing.
    
    This class extends OutputWriter to handle multiple batches efficiently,
    managing file handles across batch boundaries and providing progress tracking.
    """
    
    def __init__(self, config: ProcessingConfig, result: ProcessingResult) -> None:
        """Initialize BatchOutputWriter.
        
        Args:
            config: Processing configuration
            result: Processing result to track output files
        """
        super().__init__(config, result)
        self._batch_count = 0
        self._total_rows_written = 0
    
    def start_batch_processing(self) -> None:
        """Start batch processing mode by opening files."""
        self.__enter__()
    
    def finish_batch_processing(self) -> None:
        """Finish batch processing mode by closing files."""
        self.__exit__(None, None, None)
        logger.info(f"Batch processing completed: {self._batch_count} batches, "
                   f"{self._total_rows_written} total rows written")
    
    def write_batch(self, gt_df: pd.DataFrame, seq_df: pd.DataFrame) -> None:
        """Write a batch with progress tracking.
        
        Args:
            gt_df: DataFrame containing genotype data
            seq_df: DataFrame containing sequence data
        """
        super().write_batch(gt_df, seq_df)
        
        self._batch_count += 1
        self._total_rows_written += len(gt_df)
        
        if self._batch_count % 10 == 0:  # Log progress every 10 batches
            logger.info(f"Processed {self._batch_count} batches, "
                       f"{self._total_rows_written} rows written")
    
    @property
    def batch_count(self) -> int:
        """Get number of batches processed."""
        return self._batch_count
    
    @property
    def total_rows_written(self) -> int:
        """Get total number of rows written."""
        return self._total_rows_written