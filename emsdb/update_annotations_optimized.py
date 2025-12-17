#!/usr/bin/env python3
"""
优化的注释更新工具 - 内存高效的基因组注释处理系统

这个工具比较注释文件和目标文件，识别注释不一致的记录，
并输出更新后的结果。使用内存优化技术处理大型文件。

作者: AI Assistant
版本: 2.0 (优化版)
"""

import sys
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console

from config import ProcessingConfig
from logger import AnnotationLogger
from validation import FileValidator, ValidationError
from processor import MemoryOptimizedProcessor
from ultra_processor import UltraMemoryProcessor
from fast_processor import FastMemoryProcessor

console = Console()

# Create the typer app
app = typer.Typer(
    name="update-annotations",
    help="优化的注释更新工具 - 处理大型基因组注释文件的内存高效工具",
    add_completion=False
)


@app.command()
def main(
    annotation_file: Path = typer.Argument(
        ...,
        help="注释文件路径 (.txt 格式，制表符分隔)",
        metavar="ANNOTATION_FILE"
    ),
    target_file: Path = typer.Argument(
        ..., 
        help="目标文件路径 (.tsv 或 .tsv.gz 格式)",
        metavar="TARGET_FILE"
    ),
    output_file: Path = typer.Argument(
        ...,
        help="输出文件路径 (.tsv 格式)",
        metavar="OUTPUT_FILE"
    ),
    chunk_size: int = typer.Option(
        10000,
        "--chunk-size", "-c",
        help="处理块大小，用于内存优化",
        min=1000,
        max=1000000
    ),
    log_level: str = typer.Option(
        "INFO",
        "--log-level", "-l",
        help="日志级别 (DEBUG, INFO, WARNING, ERROR)",
        case_sensitive=False
    ),
    memory_limit: Optional[float] = typer.Option(
        None,
        "--memory-limit", "-m",
        help="内存限制 (GB)，超过此限制将调整处理策略",
        min=0.1
    ),
    log_file: Optional[Path] = typer.Option(
        None,
        "--log-file",
        help="可选的日志文件路径"
    ),
    ultra_mode: bool = typer.Option(
        False,
        "--ultra-memory",
        help="启用超级内存优化模式（分批处理，适用于内存限制严格的环境）"
    )
):
    """
    处理基因组注释文件，识别并更新不一致的记录。
    
    此工具比较注释文件和目标文件，找出注释不一致的记录，
    并输出更新后的结果。使用内存优化技术处理大型文件。
    """
    
    # Validate log level
    valid_log_levels = ["DEBUG", "INFO", "WARNING", "ERROR"]
    log_level = log_level.upper()
    if log_level not in valid_log_levels:
        console.print(f"[red]错误: 无效的日志级别 '{log_level}'. 有效选项: {', '.join(valid_log_levels)}[/red]")
        raise typer.Exit(1)
    
    # Validate input files exist
    for file_path, file_type in [(annotation_file, "注释"), (target_file, "目标")]:
        if not file_path.exists():
            console.print(f"[red]错误: {file_type}文件不存在: {file_path}[/red]")
            raise typer.Exit(1)
        
        if not file_path.is_file():
            console.print(f"[red]错误: {file_path} 不是一个文件[/red]")
            raise typer.Exit(1)
    
    # Validate output directory
    output_dir = output_file.parent
    if not output_dir.exists():
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            console.print(f"[red]错误: 无法创建输出目录 {output_dir}: {str(e)}[/red]")
            raise typer.Exit(1)
    
    # Create configuration
    config = ProcessingConfig(
        annotation_file=annotation_file,
        target_file=target_file,
        output_file=output_file,
        chunk_size=chunk_size,
        log_level=log_level,
        memory_limit_gb=memory_limit
    )
    
    # Run the processing pipeline
    run_processing_pipeline(config)


def run_processing_pipeline(config: ProcessingConfig):
    """Run the complete processing pipeline with the given configuration."""
    
    # Initialize logging system
    annotation_logger = AnnotationLogger(config)
    
    try:
        # Log startup information
        annotation_logger.log_startup()
        annotation_logger.log_system_info()
        
        # Step 1: Validate input files
        annotation_logger.log_step_start("文件验证", 1, 3)
        validator = FileValidator(config)
        file_metadata = validator.validate_all_files()
        annotation_logger.log_step_complete("文件验证")
        
        # Log file information
        for file_type, metadata in file_metadata.items():
            annotation_logger.log_file_info(metadata.path, file_type)
        
        # Step 2: Initialize processor (choose based on memory requirements)
        annotation_logger.log_step_start("初始化处理器", 2, 3)
        
        # Choose processor based on memory requirements and file size
        total_size_gb = sum(m.size_bytes for m in file_metadata.values()) / (1024**3)
        
        if config.memory_limit_gb and config.memory_limit_gb <= 20:
            # Ultra mode for very strict memory limits
            processor = UltraMemoryProcessor(config, annotation_logger)
            annotation_logger.log_progress("使用超级内存优化处理器（分批处理模式）")
        elif config.memory_limit_gb and config.memory_limit_gb <= 30:
            # Fast mode for 30GB limit
            processor = FastMemoryProcessor(config, annotation_logger)
            annotation_logger.log_progress("使用快速内存优化处理器（30GB优化模式）")
        else:
            # Standard mode for higher memory limits
            processor = MemoryOptimizedProcessor(config, annotation_logger)
            annotation_logger.log_progress("使用标准内存优化处理器")
        
        # Optimize chunk size based on file sizes and memory constraints
        total_size_mb = sum(m.size_bytes for m in file_metadata.values()) / (1024 * 1024)
        
        # If no memory limit is set and files are large, set a reasonable limit
        if not config.memory_limit_gb and total_size_mb > 3000:  # Files > 3GB
            config.memory_limit_gb = 30.0  # 30GB limit for better performance
            annotation_logger.log_progress(f"检测到大文件 ({total_size_mb:.1f} MB)，自动设置内存限制为 {config.memory_limit_gb} GB")
        
        # Only optimize chunk size for standard processor
        if hasattr(processor, 'optimize_chunk_size_for_memory'):
            optimized_chunk_size = processor.optimize_chunk_size_for_memory(total_size_mb)
            if optimized_chunk_size != config.chunk_size:
                config.chunk_size = optimized_chunk_size
        
        annotation_logger.log_step_complete("初始化处理器")
        
        # Step 3: Process annotations
        annotation_logger.log_step_start("处理注释数据", 3, 3)
        inconsistent_count = processor.process_annotations(file_metadata)
        annotation_logger.log_step_complete("处理注释数据")
        
        # Log completion summary
        annotation_logger.log_completion(inconsistent_count)
        
        # Display success message
        console.print(f"\n[green]✅ 处理成功完成！[/green]")
        console.print(f"[cyan]发现 {inconsistent_count:,} 条不一致记录[/cyan]")
        console.print(f"[cyan]结果已保存到: {config.output_file}[/cyan]")
        
    except ValidationError as e:
        annotation_logger.log_error(e, "文件验证")
        console.print(f"\n[red]❌ 验证失败: {str(e)}[/red]")
        console.print("[yellow]请检查输入文件格式和内容[/yellow]")
        sys.exit(1)
        
    except Exception as e:
        annotation_logger.log_error(e, "处理流程")
        console.print(f"\n[red]❌ 处理失败: {str(e)}[/red]")
        console.print("[yellow]请查看日志获取详细错误信息[/yellow]")
        sys.exit(1)


def run_with_args(
    annotation_file: str,
    target_file: str, 
    output_file: str,
    chunk_size: int = 10000,
    log_level: str = "INFO",
    memory_limit: Optional[float] = None
):
    """
    Programmatic interface for running the processor.
    
    Args:
        annotation_file: Path to annotation file
        target_file: Path to target file
        output_file: Path to output file
        chunk_size: Processing chunk size
        log_level: Logging level
        memory_limit: Memory limit in GB
    """
    config = ProcessingConfig(
        annotation_file=Path(annotation_file),
        target_file=Path(target_file),
        output_file=Path(output_file),
        chunk_size=chunk_size,
        log_level=log_level,
        memory_limit_gb=memory_limit
    )
    
    run_processing_pipeline(config)


if __name__ == "__main__":
    app()