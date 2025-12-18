"""Input file validation system for the annotation processing pipeline."""

import csv
import gzip
from pathlib import Path
from typing import List, Dict, Set, Optional, Tuple
import polars as pl
from loguru import logger

from config import ProcessingConfig, FileMetadata, ANNOTATION_COLUMNS, JOIN_COLUMNS


class ValidationError(Exception):
    """Custom exception for validation errors."""
    pass


class FileValidator:
    """Validates input files for format, schema, and accessibility."""
    
    def __init__(self, config: ProcessingConfig):
        self.config = config
    
    def validate_all_files(self) -> Dict[str, FileMetadata]:
        """Validate all input files and return metadata."""
        logger.info("开始验证输入文件...")
        
        metadata = {}
        
        try:
            # Validate annotation file
            logger.info("验证注释文件...")
            metadata['annotation'] = self.validate_annotation_file(self.config.annotation_file)
            
            # Validate target file  
            logger.info("验证目标文件...")
            metadata['target'] = self.validate_target_file(self.config.target_file)
            
            # Cross-validate files
            self.cross_validate_files(metadata['annotation'], metadata['target'])
            
            logger.success("所有文件验证通过")
            return metadata
            
        except ValidationError as e:
            logger.error(f"文件验证失败: {str(e)}")
            raise
        except Exception as e:
            logger.error(f"验证过程中发生意外错误: {str(e)}")
            raise ValidationError(f"验证失败: {str(e)}")
    
    def validate_annotation_file(self, file_path: Path) -> FileMetadata:
        """Validate annotation file format and schema."""
        if not file_path.exists():
            raise ValidationError(f"注释文件不存在: {file_path}")
        
        # Check file size
        size_bytes = file_path.stat().st_size
        if size_bytes == 0:
            raise ValidationError(f"注释文件为空: {file_path}")
        
        # Check if file is readable
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                first_line = f.readline().strip()
        except UnicodeDecodeError:
            try:
                with open(file_path, 'r', encoding='latin-1') as f:
                    first_line = f.readline().strip()
                logger.warning("注释文件使用 latin-1 编码")
            except Exception as e:
                raise ValidationError(f"无法读取注释文件 {file_path}: {str(e)}")
        except Exception as e:
            raise ValidationError(f"无法访问注释文件 {file_path}: {str(e)}")
        
        # Validate format (should be tab-separated)
        if not first_line:
            raise ValidationError(f"注释文件第一行为空: {file_path}")
        
        # Count actual columns including empty ones at the end
        # Use awk-style counting to handle trailing tabs correctly
        import subprocess
        try:
            result = subprocess.run(
                ["awk", "-F\\t", "{print NF; exit}", str(file_path)],
                capture_output=True, text=True, check=True
            )
            actual_columns = int(result.stdout.strip())
        except:
            # Fallback to Python split if awk fails
            actual_columns = len(first_line.split('\t'))
        
        if actual_columns < len(ANNOTATION_COLUMNS):
            raise ValidationError(
                f"注释文件列数不足: 期望至少 {len(ANNOTATION_COLUMNS)} 列，实际 {actual_columns} 列"
            )
        
        logger.debug(f"注释文件实际列数: {actual_columns}, 期望: {len(ANNOTATION_COLUMNS)}")
        columns = first_line.split('\t')
        # Pad columns list to match actual column count
        while len(columns) < actual_columns:
            columns.append("")
        
        # Try to read a sample with polars to validate format
        try:
            # Read all columns first, then select only the first 10
            # Create column names for all columns in the file
            all_column_names = ANNOTATION_COLUMNS.copy()
            if actual_columns > len(ANNOTATION_COLUMNS):
                # Add names for extra columns
                for i in range(len(ANNOTATION_COLUMNS), actual_columns):
                    all_column_names.append(f"extra_{i}")
            
            sample_df = pl.read_csv(
                file_path,
                separator='\t',
                has_header=False,
                n_rows=100,
                new_columns=all_column_names
            ).select(ANNOTATION_COLUMNS)  # Only keep the first 10 columns
            
            # Validate required columns can be parsed
            for col in JOIN_COLUMNS:
                if col in sample_df.columns:
                    if col == 'pos':
                        # Check if position column contains valid integers
                        try:
                            sample_df.select(pl.col(col).cast(pl.Int64))
                        except Exception:
                            raise ValidationError(f"注释文件中 '{col}' 列包含无效的数值")
        
        except Exception as e:
            if isinstance(e, ValidationError):
                raise
            raise ValidationError(f"注释文件格式错误: {str(e)}")
        
        # Estimate row count
        row_count = self._estimate_row_count(file_path)
        
        return FileMetadata(
            path=file_path,
            size_bytes=size_bytes,
            row_count=row_count,
            column_count=len(columns),
            schema={col: "String" for col in ANNOTATION_COLUMNS}
        )
    
    def validate_target_file(self, file_path: Path) -> FileMetadata:
        """Validate target file format and schema."""
        if not file_path.exists():
            raise ValidationError(f"目标文件不存在: {file_path}")
        
        # Check file size
        size_bytes = file_path.stat().st_size
        if size_bytes == 0:
            raise ValidationError(f"目标文件为空: {file_path}")
        
        # Handle compressed files
        is_compressed = file_path.suffix.lower() == '.gz'
        
        try:
            # Read header to validate format
            if is_compressed:
                with gzip.open(file_path, 'rt', encoding='utf-8') as f:
                    header_line = f.readline().strip()
            else:
                with open(file_path, 'r', encoding='utf-8') as f:
                    header_line = f.readline().strip()
        except Exception as e:
            raise ValidationError(f"无法读取目标文件 {file_path}: {str(e)}")
        
        if not header_line:
            raise ValidationError(f"目标文件头部为空: {file_path}")
        
        # Parse header
        columns = header_line.split('\t')
        
        # Check for required join columns
        missing_columns = []
        for col in JOIN_COLUMNS:
            if col not in columns:
                missing_columns.append(col)
        
        if missing_columns:
            raise ValidationError(
                f"目标文件缺少必需的列: {', '.join(missing_columns)}. "
                f"可用列: {', '.join(columns)}"
            )
        
        # Try to read a sample with polars
        try:
            sample_df = pl.read_csv(
                file_path,
                separator='\t',
                has_header=True,
                n_rows=100
            )
            
            # Validate join columns
            for col in JOIN_COLUMNS:
                if col == 'pos':
                    try:
                        sample_df.select(pl.col(col).cast(pl.Int64))
                    except Exception:
                        raise ValidationError(f"目标文件中 '{col}' 列包含无效的数值")
        
        except Exception as e:
            if isinstance(e, ValidationError):
                raise
            raise ValidationError(f"目标文件格式错误: {str(e)}")
        
        # Estimate row count
        row_count = self._estimate_row_count(file_path, is_compressed)
        
        return FileMetadata(
            path=file_path,
            size_bytes=size_bytes,
            row_count=row_count,
            column_count=len(columns),
            schema={col: str(sample_df[col].dtype) for col in columns}
        )
    
    def cross_validate_files(self, ann_metadata: FileMetadata, target_metadata: FileMetadata):
        """Cross-validate that files are compatible."""
        # Check that both files have reasonable sizes
        min_size = 1024  # 1KB minimum
        
        if ann_metadata.size_bytes < min_size:
            raise ValidationError("注释文件太小，可能不包含有效数据")
        
        if target_metadata.size_bytes < min_size:
            raise ValidationError("目标文件太小，可能不包含有效数据")
        
        # Log file information
        logger.info(f"注释文件: {ann_metadata.row_count:,} 行, {ann_metadata.size_bytes / (1024*1024):.1f} MB")
        logger.info(f"目标文件: {target_metadata.row_count:,} 行, {target_metadata.size_bytes / (1024*1024):.1f} MB")
        
        # Warn if files are very large
        large_file_threshold = 1024 * 1024 * 1024  # 1GB
        
        if ann_metadata.size_bytes > large_file_threshold:
            logger.warning(f"注释文件很大 ({ann_metadata.size_bytes / (1024*1024*1024):.1f} GB)，处理可能需要较长时间")
        
        if target_metadata.size_bytes > large_file_threshold:
            logger.warning(f"目标文件很大 ({target_metadata.size_bytes / (1024*1024*1024):.1f} GB)，处理可能需要较长时间")
    
    def _estimate_row_count(self, file_path: Path, is_compressed: bool = False) -> Optional[int]:
        """Estimate the number of rows in a file."""
        try:
            if is_compressed:
                # For compressed files, read a sample and estimate
                with gzip.open(file_path, 'rt') as f:
                    sample_lines = 0
                    sample_bytes = 0
                    for line in f:
                        sample_lines += 1
                        sample_bytes += len(line.encode('utf-8'))
                        if sample_lines >= 1000:  # Sample first 1000 lines
                            break
                
                if sample_lines > 0:
                    # Estimate total lines based on sample
                    total_bytes = file_path.stat().st_size
                    estimated_rows = int((total_bytes / sample_bytes) * sample_lines)
                    return max(1, estimated_rows - 1)  # Subtract header
            else:
                # For uncompressed files, count lines more efficiently
                with open(file_path, 'rb') as f:
                    line_count = sum(1 for _ in f)
                return max(1, line_count - 1)  # Subtract header
        
        except Exception as e:
            logger.warning(f"无法估算文件行数 {file_path}: {str(e)}")
            return None