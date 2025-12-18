"""Ultra memory-optimized data processing with batch processing."""

import gc
import tempfile
from pathlib import Path
from typing import Dict, Optional, List
import polars as pl
from loguru import logger

from config import ProcessingConfig, FileMetadata, ANNOTATION_COLUMNS, UPDATE_COLUMNS, JOIN_COLUMNS
from logger import AnnotationLogger


class UltraMemoryProcessor:
    """Ultra memory-efficient processor using batch processing."""
    
    def __init__(self, config: ProcessingConfig, annotation_logger: AnnotationLogger):
        self.config = config
        self.logger = annotation_logger
        self.temp_dir = Path(tempfile.mkdtemp(prefix="annotation_processing_"))
        self._setup_polars_config()
    
    def _setup_polars_config(self):
        """Configure Polars for ultra memory optimization."""
        # Use extremely small chunk size
        ultra_chunk_size = 500  # Very small chunks
        pl.Config.set_streaming_chunk_size(ultra_chunk_size)
        
        logger.info(f"启用超级内存优化模式，块大小: {ultra_chunk_size:,}")
        logger.info(f"临时目录: {self.temp_dir}")
    
    def process_annotations(self, file_metadata: Dict[str, FileMetadata]) -> int:
        """Process annotation files using batch processing to minimize memory."""
        progress = self.logger.create_progress_tracker(6)
        
        try:
            # Step 1: Create annotation index for faster lookups
            progress.start_step("创建注释索引")
            ann_index_file = self._create_annotation_index()
            progress.complete_step("创建注释索引")
            
            # Step 2: Get target file schema
            progress.start_step("分析目标文件结构")
            target_cols = self._get_target_schema()
            progress.complete_step("分析目标文件结构")
            
            # Step 3: Process target file in batches
            progress.start_step("分批处理目标文件")
            batch_files = self._process_target_in_batches(ann_index_file, target_cols)
            progress.complete_step("分批处理目标文件")
            
            # Step 4: Merge batch results
            progress.start_step("合并批处理结果")
            self._merge_batch_results(batch_files, target_cols)
            progress.complete_step("合并批处理结果")
            
            # Step 5: Count results
            progress.start_step("统计结果")
            inconsistent_count = self._count_final_results()
            progress.complete_step("统计结果")
            
            # Step 6: Cleanup
            progress.start_step("清理临时文件")
            self._cleanup_temp_files()
            progress.complete_step("清理临时文件")
            
            return inconsistent_count
            
        except Exception as e:
            self.logger.log_error(e, "超级内存优化处理")
            raise
        finally:
            # Cleanup temp directory
            self._cleanup_temp_files()
            gc.collect()
    
    def _create_annotation_index(self) -> Path:
        """Create a memory-efficient annotation index."""
        index_file = self.temp_dir / "annotation_index.parquet"
        
        # Read annotation file in small chunks and create index
        chunk_size = 10000  # Small chunks for indexing
        
        logger.info("正在创建注释索引...")
        
        # Process annotation file in chunks
        chunks_processed = 0
        index_chunks = []
        
        # Read in chunks using scan_csv with collect in batches
        # Step 1: Scan without new_columns
        lf_ann = pl.scan_csv(
            self.config.annotation_file,
            separator="\t",
            has_header=False,
            low_memory=True,
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
            # Only keep essential columns for indexing
            pl.col("chrom").cast(pl.Categorical),
            pl.col("pos").cast(pl.Int32),
            pl.col("refer").cast(pl.Categorical),
            pl.col("alt").cast(pl.Categorical),
            pl.col("type").cast(pl.Categorical),
            pl.col("impact").cast(pl.Categorical),
            pl.col("gene").cast(pl.Utf8),
            pl.col("exon_rank").cast(pl.Utf8),
            pl.col("cds_pos").cast(pl.Utf8),
            pl.col("protein_pos").cast(pl.Utf8),
        ])
        
        # Write to parquet for efficient access
        lf_ann.sink_parquet(index_file)
        
        logger.info(f"注释索引已创建: {index_file}")
        return index_file
    
    def _get_target_schema(self) -> List[str]:
        """Get target file schema with minimal memory usage."""
        sample_df = pl.read_csv(
            self.config.target_file,
            separator="\t",
            has_header=True,
            n_rows=1
        )
        target_cols = sample_df.columns
        del sample_df
        gc.collect()
        return target_cols
    
    def _process_target_in_batches(self, ann_index_file: Path, target_cols: List[str]) -> List[Path]:
        """Process target file in small batches to control memory usage."""
        batch_size = 5000  # Very small batches
        batch_files = []
        batch_num = 0
        
        logger.info(f"开始分批处理，批大小: {batch_size:,}")
        
        # Load annotation index once
        lf_ann_index = pl.scan_parquet(ann_index_file)
        
        # Process target file in batches
        lf_target = pl.scan_csv(
            self.config.target_file,
            separator="\t",
            has_header=True,
            low_memory=True,
        ).with_columns([
            pl.col("chrom").cast(pl.Categorical),
            pl.col("pos").cast(pl.Int32),
            pl.col("refer").cast(pl.Categorical),
            pl.col("alt").cast(pl.Categorical),
        ])
        
        # Get total rows for progress tracking
        try:
            total_rows = lf_target.select(pl.len()).collect().item()
            logger.info(f"目标文件总行数: {total_rows:,}")
        except:
            total_rows = None
        
        # Process in batches
        offset = 0
        while True:
            try:
                # Read a batch
                batch_lf = lf_target.slice(offset, batch_size)
                
                # Check if batch is empty
                batch_df = batch_lf.collect()
                if len(batch_df) == 0:
                    break
                
                logger.debug(f"处理批次 {batch_num + 1}，行数: {len(batch_df):,}")
                
                # Convert back to lazy for processing
                batch_lf = batch_df.lazy()
                del batch_df
                
                # Join with annotation index
                joined_lf = batch_lf.join(
                    lf_ann_index,
                    on=JOIN_COLUMNS,
                    how="inner",
                    suffix="_new"
                )
                
                # Filter for inconsistent records
                # Compare in priority order: cds_pos, protein_pos, type
                filter_conditions = []
                comparison_cols = ["cds_pos", "protein_pos", "type"]  # Priority order
                for col in comparison_cols:
                    if col in target_cols:
                        # Must cast to Utf8 because Categorical columns cannot be coalesced with string literals
                        condition = (
                            pl.coalesce([pl.col(col).cast(pl.Utf8), pl.lit("")]) != 
                            pl.coalesce([pl.col(f"{col}_new").cast(pl.Utf8), pl.lit("")])
                        )
                        filter_conditions.append(condition)
                
                if filter_conditions:
                    combined_filter = filter_conditions[0]
                    for condition in filter_conditions[1:]:
                        combined_filter = combined_filter | condition
                    
                    filtered_lf = joined_lf.filter(combined_filter)
                    
                    # Update annotation columns
                    select_exprs = []
                    for col in target_cols:
                        if col in UPDATE_COLUMNS:
                            new_col = f"{col}_new"
                            select_exprs.append(
                                pl.when(pl.col(new_col).is_not_null() & (pl.col(new_col) != ""))
                                .then(pl.col(new_col))
                                .otherwise(pl.col(col))
                                .alias(col)
                            )
                        else:
                            select_exprs.append(pl.col(col))
                    
                    final_lf = filtered_lf.select(select_exprs)
                    
                    # Save batch result if not empty
                    batch_result = final_lf.collect()
                    if len(batch_result) > 0:
                        batch_file = self.temp_dir / f"batch_{batch_num:06d}.parquet"
                        batch_result.write_parquet(batch_file)
                        batch_files.append(batch_file)
                        logger.debug(f"批次 {batch_num + 1} 保存了 {len(batch_result):,} 条不一致记录")
                    
                    del batch_result, final_lf, filtered_lf
                
                del joined_lf, batch_lf
                gc.collect()
                
                # Monitor memory
                self.logger.monitor_memory()
                
                batch_num += 1
                offset += batch_size
                
                # Progress update
                if total_rows and batch_num % 10 == 0:
                    progress_pct = min(100, (offset / total_rows) * 100)
                    logger.info(f"已处理 {batch_num} 个批次 ({progress_pct:.1f}%)")
                
            except Exception as e:
                logger.error(f"处理批次 {batch_num} 时出错: {str(e)}")
                break
        
        logger.info(f"分批处理完成，共处理 {batch_num} 个批次，生成 {len(batch_files)} 个结果文件")
        return batch_files
    
    def _merge_batch_results(self, batch_files: List[Path], target_cols: List[str]):
        """Merge batch results into final output file."""
        if not batch_files:
            # Create empty output file
            empty_df = pl.DataFrame({col: [] for col in target_cols})
            empty_df.write_csv(self.config.output_file, separator="\t")
            return
        
        logger.info(f"合并 {len(batch_files)} 个批次结果文件...")
        
        # Read and concatenate all batch files
        batch_dfs = []
        for i, batch_file in enumerate(batch_files):
            try:
                batch_df = pl.read_parquet(batch_file)
                batch_dfs.append(batch_df)
                
                # Merge in chunks to control memory
                if len(batch_dfs) >= 10:  # Merge every 10 files
                    merged_df = pl.concat(batch_dfs)
                    batch_dfs = [merged_df]
                    gc.collect()
                    
                logger.debug(f"已读取批次文件 {i + 1}/{len(batch_files)}")
                
            except Exception as e:
                logger.warning(f"读取批次文件 {batch_file} 失败: {str(e)}")
        
        # Final merge
        if batch_dfs:
            final_df = pl.concat(batch_dfs)
            
            # Write to output file
            final_df.write_csv(self.config.output_file, separator="\t")
            
            logger.info(f"合并完成，输出 {len(final_df):,} 条记录")
            del final_df
        
        gc.collect()
    
    def _count_final_results(self) -> int:
        """Count records in final output file."""
        try:
            import subprocess
            result = subprocess.run(
                ["wc", "-l", str(self.config.output_file)],
                capture_output=True, text=True, check=True
            )
            total_lines = int(result.stdout.split()[0])
            return max(0, total_lines - 1)  # Subtract header
        except:
            try:
                df = pl.read_csv(self.config.output_file, separator="\t", n_rows=0)
                lf = pl.scan_csv(self.config.output_file, separator="\t", has_header=True)
                return lf.select(pl.len()).collect().item()
            except:
                return 0
    
    def _cleanup_temp_files(self):
        """Clean up temporary files."""
        try:
            import shutil
            if self.temp_dir.exists():
                shutil.rmtree(self.temp_dir)
                logger.debug(f"已清理临时目录: {self.temp_dir}")
        except Exception as e:
            logger.warning(f"清理临时文件失败: {str(e)}")