"""Memory-optimized data processing using Polars lazy evaluation."""

import gc
from pathlib import Path
from typing import Dict, Optional
import polars as pl
from loguru import logger

from config import ProcessingConfig, FileMetadata, ANNOTATION_COLUMNS, UPDATE_COLUMNS, JOIN_COLUMNS
from logger import AnnotationLogger


class MemoryOptimizedProcessor:
    """Memory-efficient processor for large annotation files."""
    
    def __init__(self, config: ProcessingConfig, annotation_logger: AnnotationLogger):
        self.config = config
        self.logger = annotation_logger
        self._setup_polars_config()
    
    def _setup_polars_config(self):
        """Configure Polars for aggressive memory optimization."""
        # Use much smaller chunk size for extreme memory efficiency
        aggressive_chunk_size = min(self.config.chunk_size, 5000)  # Max 5000 rows per chunk
        pl.Config.set_streaming_chunk_size(aggressive_chunk_size)
        
        # Enable streaming mode (if available in this Polars version)
        try:
            pl.Config.set_streaming_engine(True)
        except AttributeError:
            # Older Polars version doesn't have this method
            logger.debug("Polars版本不支持set_streaming_engine，跳过此设置")
        
        # Configure memory settings if limit is specified
        if self.config.memory_limit_gb:
            # Convert GB to bytes for internal calculations
            memory_limit_bytes = int(self.config.memory_limit_gb * 1024 * 1024 * 1024)
            logger.debug(f"配置内存限制: {self.config.memory_limit_gb:.1f} GB")
            # Use even smaller chunks if memory limit is set
            if self.config.memory_limit_gb <= 20:
                aggressive_chunk_size = min(aggressive_chunk_size, 2000)
                pl.Config.set_streaming_chunk_size(aggressive_chunk_size)
                logger.info(f"启用极限内存模式，块大小: {aggressive_chunk_size:,}")
        
        logger.debug(f"Polars流式块大小: {aggressive_chunk_size:,}")
    
    def process_annotations(self, file_metadata: Dict[str, FileMetadata]) -> int:
        """Process annotation files with aggressive memory optimization."""
        progress = self.logger.create_progress_tracker(4)
        
        try:
            # Monitor initial memory
            initial_memory = self.logger.monitor_memory()
            
            # Step 1: Load annotation file lazily
            progress.start_step("加载注释文件")
            lf_ann = self._load_annotation_file_lazy()
            progress.complete_step("加载注释文件")
            
            # Force garbage collection after each step
            gc.collect()
            self.logger.monitor_memory()
            
            # Step 2: Load target file lazily
            progress.start_step("加载目标文件")
            lf_target, target_cols = self._load_target_file_lazy()
            progress.complete_step("加载目标文件")
            
            # Force garbage collection
            gc.collect()
            self.logger.monitor_memory()
            
            # Step 3: Perform memory-efficient join and filtering
            progress.start_step("比较并筛选不一致记录")
            lf_filtered = self._perform_memory_efficient_join(lf_ann, lf_target, target_cols)
            progress.complete_step("比较并筛选不一致记录")
            
            # Clear references to large objects
            del lf_ann, lf_target
            gc.collect()
            self.logger.monitor_memory()
            
            # Step 4: Stream results to output file
            progress.start_step("写入结果文件")
            inconsistent_count = self._stream_results_to_file(lf_filtered)
            progress.complete_step("写入结果文件")
            
            # Final cleanup
            del lf_filtered
            gc.collect()
            final_memory = self.logger.monitor_memory()
            
            logger.info(f"内存使用变化: {initial_memory:.1f} MB -> {final_memory:.1f} MB")
            
            return inconsistent_count
            
        except Exception as e:
            self.logger.log_error(e, "数据处理")
            raise
        finally:
            # Aggressive garbage collection to free memory
            gc.collect()
            gc.collect()  # Run twice for better cleanup
    
    def _load_annotation_file_lazy(self) -> pl.LazyFrame:
        """Load annotation file using aggressive memory optimization."""
        try:
            # Use very small infer_schema_length to reduce memory usage
            small_infer_length = min(1000, self.config.chunk_size // 10)
            
            # Step 1: Scan without new_columns to let Polars handle parsing naturally
            lf_ann = pl.scan_csv(
                self.config.annotation_file,
                separator="\t",
                has_header=False,
                infer_schema_length=small_infer_length,  # Much smaller schema inference
                low_memory=True,  # Enable low memory mode
                rechunk=False,  # Disable rechunking to save memory
            )
            
            # Step 2: Dynamically rename columns based on position
            current_cols = lf_ann.collect_schema().names()
            
            rename_map = {}
            for i, col_name in enumerate(current_cols):
                if i < len(ANNOTATION_COLUMNS):
                    rename_map[col_name] = ANNOTATION_COLUMNS[i]
                else:
                    rename_map[col_name] = f"extra_{i}"
            
            lf_ann = lf_ann.rename(rename_map)
            
            # Step 3: Select and cast
            lf_ann = lf_ann.select([
                # Use more memory-efficient data types
                pl.col("chrom").cast(pl.Categorical),  # Categorical for repeated values
                pl.col("pos").cast(pl.Int32),  # Use Int32 instead of Int64 if positions fit
                pl.col("refer").cast(pl.Categorical),  # Categorical for DNA bases
                pl.col("alt").cast(pl.Categorical),   # Categorical for DNA bases
                pl.col("type").cast(pl.Categorical),
                pl.col("impact").cast(pl.Categorical),
                pl.col("gene").cast(pl.Utf8),
                pl.col("exon_rank").cast(pl.Utf8),
                pl.col("cds_pos").cast(pl.Utf8),
                pl.col("protein_pos").cast(pl.Utf8),
            ])
            
            logger.debug("注释文件已配置为极限内存优化模式")
            return lf_ann
            
        except Exception as e:
            raise Exception(f"加载注释文件失败: {str(e)}")
    
    def _load_target_file_lazy(self) -> tuple[pl.LazyFrame, list]:
        """Load target file using aggressive memory optimization."""
        try:
            # First, get the schema to preserve column order (minimal read)
            sample_df = pl.read_csv(
                self.config.target_file,
                separator="\t",
                has_header=True,
                n_rows=1
            )
            target_cols = sample_df.columns
            del sample_df  # Immediately free memory
            gc.collect()
            
            # Use very small infer_schema_length to reduce memory usage
            small_infer_length = min(1000, self.config.chunk_size // 10)
            
            # Load lazily with aggressive memory optimization
            lf_target = pl.scan_csv(
                self.config.target_file,
                separator="\t",
                has_header=True,
                infer_schema_length=small_infer_length,  # Much smaller schema inference
                low_memory=True,
                rechunk=False,  # Disable rechunking to save memory
            ).with_columns([
                # Use memory-efficient data types for join columns
                pl.col("chrom").cast(pl.Categorical),
                pl.col("pos").cast(pl.Int32),  # Use Int32 instead of Int64
                pl.col("refer").cast(pl.Categorical),  # Categorical for DNA bases
                pl.col("alt").cast(pl.Categorical),    # Categorical for DNA bases
            ])
            
            logger.debug("目标文件已配置为极限内存优化模式")
            return lf_target, target_cols
            
        except Exception as e:
            raise Exception(f"加载目标文件失败: {str(e)}")
    
    def _perform_memory_efficient_join(self, lf_ann: pl.LazyFrame, lf_target: pl.LazyFrame, target_cols: list) -> pl.LazyFrame:
        """Perform ultra memory-efficient join and filtering with streaming."""
        try:
            # Use streaming join with minimal memory footprint
            lf_joined = lf_target.join(
                lf_ann,
                on=JOIN_COLUMNS,
                how="inner",
                suffix="_new",
                join_nulls=False,  # Skip null joins to save memory
            )
            
            # Filter for inconsistent records using lazy evaluation
            # Compare in priority order: cds_pos, protein_pos, type
            filter_conditions = []
            comparison_cols = ["cds_pos", "protein_pos", "type"]  # Priority order
            for col in comparison_cols:
                if col in target_cols:  # Only filter if column exists in target
                    # Use coalesce to handle nulls efficiently
                    # Must cast to Utf8 because Categorical columns cannot be coalesced with string literals
                    condition = (
                        pl.coalesce([pl.col(col).cast(pl.Utf8), pl.lit("")]) != 
                        pl.coalesce([pl.col(f"{col}_new").cast(pl.Utf8), pl.lit("")])
                    )
                    filter_conditions.append(condition)
            
            if not filter_conditions:
                raise Exception("目标文件中未找到可比较的注释列")
            
            # Combine conditions with OR
            combined_filter = filter_conditions[0]
            for condition in filter_conditions[1:]:
                combined_filter = combined_filter | condition
            
            lf_filtered = lf_joined.filter(combined_filter)
            
            # Update annotation columns while preserving original column order
            select_exprs = []
            for col in target_cols:
                if col in UPDATE_COLUMNS:
                    # Use new annotation if available, with null handling
                    new_col = f"{col}_new"
                    select_exprs.append(
                        pl.when(pl.col(new_col).is_not_null() & (pl.col(new_col) != ""))
                        .then(pl.col(new_col))
                        .otherwise(pl.col(col))
                        .alias(col)
                    )
                else:
                    # Keep original column
                    select_exprs.append(pl.col(col))
            
            lf_final = lf_filtered.select(select_exprs)
            
            logger.debug("连接和筛选操作已配置为极限内存优化模式")
            return lf_final
            
        except Exception as e:
            raise Exception(f"数据连接和筛选失败: {str(e)}")
    
    def _stream_results_to_file(self, lf_final: pl.LazyFrame) -> int:
        """Stream results to output file with memory monitoring."""
        try:
            # Monitor memory before streaming
            initial_memory = self.logger.monitor_memory()
            
            # Use sink_csv for streaming write (memory efficient)
            lf_final.sink_csv(
                self.config.output_file,
                separator="\t",
                include_header=True
            )
            
            # Monitor memory after streaming
            final_memory = self.logger.monitor_memory()
            memory_diff = final_memory - initial_memory
            
            if memory_diff > 0:
                logger.debug(f"写入过程内存增长: {memory_diff:.1f} MB")
            
            # Count results efficiently
            inconsistent_count = self._count_output_records()
            
            return inconsistent_count
            
        except Exception as e:
            raise Exception(f"写入结果文件失败: {str(e)}")
    
    def _count_output_records(self) -> int:
        """Efficiently count records in output file."""
        try:
            # Try using wc command first (fastest for large files)
            import subprocess
            
            result = subprocess.run(
                ["wc", "-l", str(self.config.output_file)],
                capture_output=True,
                text=True,
                check=True
            )
            total_lines = int(result.stdout.split()[0])
            return max(0, total_lines - 1)  # Subtract header
            
        except Exception:
            # Fallback to polars if wc fails
            try:
                df = pl.read_csv(
                    self.config.output_file,
                    separator="\t",
                    has_header=True,
                    n_rows=0  # Only read metadata
                )
                # Use lazy frame to count efficiently
                lf = pl.scan_csv(
                    self.config.output_file,
                    separator="\t",
                    has_header=True
                )
                return lf.select(pl.len()).collect().item()
                
            except Exception as e:
                logger.warning(f"无法统计输出记录数: {str(e)}")
                return 0
    
    def optimize_chunk_size_for_memory(self, file_size_mb: float) -> int:
        """Aggressively optimize chunk size for extreme memory efficiency."""
        # Always use aggressive optimization for large files
        if file_size_mb > 1000:  # Files larger than 1GB
            # Use very small chunks for large files
            base_chunk_size = 1000
        else:
            base_chunk_size = self.config.chunk_size
        
        if self.config.memory_limit_gb:
            memory_limit_mb = self.config.memory_limit_gb * 1024
            
            # Use extremely conservative approach: aim to use max 20% of memory limit
            target_memory_mb = memory_limit_mb * 0.2
            
            # More conservative estimate for memory per row (includes join overhead)
            estimated_bytes_per_row = 400  # Account for join operations
            max_rows_in_memory = int(target_memory_mb * 1024 * 1024 / estimated_bytes_per_row)
            
            # For 30GB limit, use reasonable chunks
            if self.config.memory_limit_gb <= 30:
                max_rows_in_memory = min(max_rows_in_memory, 5000)
                logger.info(f"启用{self.config.memory_limit_gb}GB内存限制模式")
            
            # Ensure chunk size is reasonable but small
            optimized_chunk_size = max(500, min(max_rows_in_memory, base_chunk_size))
        else:
            # Even without explicit limit, use smaller chunks for large files
            optimized_chunk_size = min(base_chunk_size, 5000)
        
        if optimized_chunk_size != self.config.chunk_size:
            logger.info(f"根据内存限制调整块大小: {self.config.chunk_size:,} -> {optimized_chunk_size:,}")
            
        return optimized_chunk_size