"""Logging utilities using loguru for comprehensive progress tracking."""

import sys
import psutil
import time
from pathlib import Path
from typing import Optional, Dict, Any
from loguru import logger
from datetime import datetime

from config import ProcessingConfig, ProcessingStats


class AnnotationLogger:
    """Centralized logging for the annotation processing system."""
    
    def __init__(self, config: ProcessingConfig, log_file: Optional[Path] = None):
        """Initialize logger with configuration."""
        self.config = config
        self.stats = ProcessingStats(start_time=datetime.now())
        
        # Remove default logger
        logger.remove()
        
        # Add console logger with appropriate level
        logger.add(
            sys.stderr,
            level=config.log_level,
            format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
            colorize=True
        )
        
        # Add file logger if specified
        if log_file:
            logger.add(
                log_file,
                level="DEBUG",
                format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {name}:{function}:{line} - {message}",
                rotation="10 MB"
            )
    
    def log_startup(self):
        """Log system startup information."""
        logger.info("=== 开始处理注释不一致的记录 ===")
        logger.info(f"注释文件: {self.config.annotation_file}")
        logger.info(f"目标文件: {self.config.target_file}")
        logger.info(f"输出文件: {self.config.output_file}")
        logger.info(f"块大小: {self.config.chunk_size:,}")
        if self.config.memory_limit_gb:
            logger.info(f"内存限制: {self.config.memory_limit_gb:.1f} GB")
    
    def log_step_start(self, step_name: str, step_number: int, total_steps: int):
        """Log the start of a processing step."""
        logger.info(f"Step {step_number}/{total_steps}: {step_name}...")
    
    def log_step_complete(self, step_name: str):
        """Log the completion of a processing step."""
        logger.success(f"  ✓ {step_name}已完成")
    
    def log_progress(self, message: str, records_processed: Optional[int] = None):
        """Log progress updates during processing."""
        if records_processed:
            self.stats.records_processed = records_processed
            logger.info(f"{message} (已处理: {records_processed:,} 条记录)")
        else:
            logger.info(message)
    
    def log_memory_usage(self, memory_mb: float):
        """Log current memory usage."""
        if memory_mb > self.stats.memory_peak_mb:
            self.stats.memory_peak_mb = memory_mb
        logger.debug(f"内存使用: {memory_mb:.1f} MB (峰值: {self.stats.memory_peak_mb:.1f} MB)")
    
    def log_error(self, error: Exception, context: str = ""):
        """Log error with context information."""
        if context:
            logger.error(f"错误在 {context}: {str(error)}")
        else:
            logger.error(f"处理错误: {str(error)}")
        logger.exception("详细错误信息:")
    
    def log_completion(self, inconsistent_count: int):
        """Log processing completion with summary statistics."""
        self.stats.end_time = datetime.now()
        self.stats.inconsistent_records = inconsistent_count
        self.stats.processing_time_seconds = (
            self.stats.end_time - self.stats.start_time
        ).total_seconds()
        
        logger.info("=" * 50)
        logger.success("处理完成！")
        logger.info(f"  不一致记录数: {inconsistent_count:,}")
        logger.info(f"  处理时间: {self.stats.processing_time_seconds:.2f} 秒")
        logger.info(f"  峰值内存: {self.stats.memory_peak_mb:.1f} MB")
        logger.info(f"  输出文件: {self.config.output_file}")
        logger.info("=" * 50)
    
    def get_stats(self) -> ProcessingStats:
        """Get current processing statistics."""
        return self.stats
    
    def log_system_info(self):
        """Log system information at startup."""
        import polars as pl
        
        # System information
        memory_gb = psutil.virtual_memory().total / (1024**3)
        cpu_count = psutil.cpu_count()
        
        logger.info(f"系统信息:")
        logger.info(f"  CPU 核心数: {cpu_count}")
        logger.info(f"  总内存: {memory_gb:.1f} GB")
        logger.info(f"  Polars 版本: {pl.__version__}")
        
        # Check available memory
        available_gb = psutil.virtual_memory().available / (1024**3)
        logger.info(f"  可用内存: {available_gb:.1f} GB")
        
        if self.config.memory_limit_gb and self.config.memory_limit_gb > available_gb:
            logger.warning(f"警告: 设置的内存限制 ({self.config.memory_limit_gb:.1f} GB) 超过可用内存")
    
    def monitor_memory(self) -> float:
        """Monitor current memory usage and return usage in MB."""
        process = psutil.Process()
        memory_info = process.memory_info()
        memory_mb = memory_info.rss / (1024 * 1024)
        
        self.log_memory_usage(memory_mb)
        
        # Check if we're approaching memory limit
        if self.config.memory_limit_gb:
            limit_mb = self.config.memory_limit_gb * 1024
            if memory_mb > limit_mb * 0.8:  # 80% of limit
                logger.warning(f"内存使用接近限制: {memory_mb:.1f} MB / {limit_mb:.1f} MB")
        
        return memory_mb
    
    def log_file_info(self, file_path: Path, file_type: str):
        """Log information about a file."""
        try:
            size_mb = file_path.stat().st_size / (1024 * 1024)
            logger.info(f"{file_type}文件大小: {size_mb:.1f} MB")
        except Exception as e:
            logger.warning(f"无法获取 {file_type}文件大小: {str(e)}")
    
    def create_progress_tracker(self, total_steps: int):
        """Create a progress tracking context."""
        return ProgressTracker(self, total_steps)


class ProgressTracker:
    """Context manager for tracking progress through multiple steps."""
    
    def __init__(self, annotation_logger: AnnotationLogger, total_steps: int):
        self.logger = annotation_logger
        self.total_steps = total_steps
        self.current_step = 0
        self.step_start_time = None
    
    def start_step(self, step_name: str):
        """Start a new processing step."""
        self.current_step += 1
        self.step_start_time = time.time()
        self.logger.log_step_start(step_name, self.current_step, self.total_steps)
        self.logger.monitor_memory()
    
    def complete_step(self, step_name: str):
        """Complete the current processing step."""
        if self.step_start_time:
            duration = time.time() - self.step_start_time
            logger.debug(f"步骤耗时: {duration:.2f} 秒")
        
        self.logger.log_step_complete(step_name)
        self.logger.monitor_memory()
    
    def update_progress(self, message: str, records_processed: Optional[int] = None):
        """Update progress within a step."""
        self.logger.log_progress(message, records_processed)