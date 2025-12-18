"""Configuration models and constants for the annotation processing system."""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Dict
from datetime import datetime


@dataclass
class ProcessingConfig:
    """Configuration for the annotation processing pipeline."""
    annotation_file: Path
    target_file: Path  
    output_file: Path
    chunk_size: int = 10000
    log_level: str = "INFO"
    memory_limit_gb: Optional[float] = None


@dataclass  
class ProcessingStats:
    """Statistics collected during processing."""
    start_time: datetime
    end_time: Optional[datetime] = None
    records_processed: int = 0
    inconsistent_records: int = 0
    memory_peak_mb: float = 0.0
    processing_time_seconds: float = 0.0


@dataclass
class FileMetadata:
    """Metadata about input/output files."""
    path: Path
    size_bytes: int
    row_count: Optional[int] = None
    column_count: int = 0
    schema: Dict[str, str] = None

    def __post_init__(self):
        if self.schema is None:
            self.schema = {}


# Default column names for annotation files (based on header.txt)
ANNOTATION_COLUMNS = [
    "chrom",
    "pos", 
    "refer",
    "alt",
    "type",
    "impact",
    "gene",
    "exon_rank",
    "cds_pos",
    "protein_pos",
]

# Columns that should be updated from annotations
UPDATE_COLUMNS = ["type", "impact", "gene", "exon_rank", "cds_pos", "protein_pos"]

# Join columns for matching records
JOIN_COLUMNS = ["chrom", "pos", "refer", "alt"]