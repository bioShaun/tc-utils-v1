#!/usr/bin/env python3
"""
VCF Processor - Standalone Script
独立的VCF处理脚本，直接使用 vcf_processor.py

这是 `python -m chip.vcf_processor.cli process` 命令的独立版本。
直接导入并使用 vcf_processor.py 中的VCFProcessor类，确保结果一致。

使用方法:
    python vcf_processor_standalone.py input.vcf output [选项]
    python vcf_processor_standalone.py input.vcf output --targets targets.txt [选项]

等价的模块化命令:
    python -m chip.vcf_processor.cli process input.vcf output [选项]
    python -m chip.vcf_processor.cli process input.vcf output --targets targets.txt [选项]

示例:
    python vcf_processor_standalone.py genotypes.vcf results
    python vcf_processor_standalone.py genotypes.vcf results --targets targets.txt
    python vcf_processor_standalone.py genotypes.vcf results --targets targets.txt --batch-size 20000 --compress

功能特性:
    - 直接使用高性能的 VCFProcessor 类，确保与模块化版本结果一致
    - 支持所有 VCF 格式（.vcf 和 .vcf.gz）
    - 可选的目标变异过滤（如果不指定则处理所有变异）
    - 批处理模式支持大文件处理
    - 完整的基因型转换和序列输出
    - 压缩输出支持
    - 详细的处理统计和进度显示

输出文件:
    output.genotype_codes.tsv: VCF格式基因型表格（CHROM, POS, REF, ALT, sample1, sample2, ...）包含 0/0, 0/1, 1/1, ./. 等编码
    output.genotype_bases.tsv: 碱基序列基因型表格（与编码表格格式相同）包含 AA, AT, TT, NN 等碱基序列

注意:
    此脚本需要相应的依赖库（cyvcf2, pandas, loguru等）。
    如果不指定目标文件，将处理VCF文件中的所有变异。
"""

import sys
import time
from pathlib import Path
from typing import Dict, Optional

import typer
from typing_extensions import Annotated

# 添加项目根目录到Python路径以支持绝对导入
current_dir = Path(__file__).parent
project_root = current_dir.parent.parent  # 从 chip/vcf_processor 到项目根目录
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))


def process_vcf(vcf_file: str, output_prefix: str, target_file: str = None,
                miss_fmt: str = "./.", gt_sep: str = "", batch_size: int = 10000,
                compress_output: bool = False, verbose: bool = False, quiet: bool = False) -> Dict:
    """使用VCFProcessor处理VCF文件
    
    Args:
        vcf_file: VCF文件路径
        output_prefix: 输出文件前缀
        target_file: 目标变异ID文件路径 (可选)
        miss_fmt: 缺失基因型格式
        gt_sep: 基因型分隔符
        batch_size: 批处理大小
        compress_output: 是否压缩输出
        verbose: 是否详细输出
        quiet: 是否静默模式
        
    Returns:
        处理结果字典
    """
    # 导入VCF处理器相关模块（在函数内部导入以确保日志已配置）
    from chip.vcf_processor.config import ProcessingConfig
    from chip.vcf_processor.vcf_processor import VCFProcessor
    
    # 创建配置
    config = ProcessingConfig(
        vcf_file=Path(vcf_file),
        output_file=Path(output_prefix),
        target_id_file=Path(target_file) if target_file else None,
        miss_fmt=miss_fmt,
        gt_sep=gt_sep,
        batch_size=batch_size,
        compress_output=compress_output,
        verbose=verbose,
        quiet=quiet,  # 正确传递静默模式参数
        threads=1  # 独立脚本使用单线程
    )
    
    # 执行处理
    start_time = time.time()
    processor = VCFProcessor(config)
    result = processor.process()
    elapsed_time = time.time() - start_time
    
    # 转换为兼容格式
    return {
        'processed_variants': result.processed_variants,
        'total_variants': result.total_variants,
        'elapsed_time': elapsed_time,
        'processing_rate': result.processed_variants / elapsed_time if elapsed_time > 0 else 0,
        'errors': result.errors,
        'success': not result.has_errors
    }


def main(
    vcf_file: Annotated[Path, typer.Argument(help="输入VCF文件路径")],
    output_prefix: Annotated[Path, typer.Argument(help="输出文件前缀")],
    targets: Annotated[Optional[Path], typer.Option("--targets", help="目标变异ID文件路径 (可选，如果不指定则处理所有变异)")] = None,
    miss_fmt: Annotated[str, typer.Option("--miss-fmt", help="缺失基因型格式")] = "./.",
    gt_sep: Annotated[str, typer.Option("--gt-sep", help="基因型分隔符")] = "",
    batch_size: Annotated[int, typer.Option("--batch-size", help="批处理大小")] = 10000,
    compress: Annotated[bool, typer.Option("--compress", help="压缩输出文件")] = False,
    verbose: Annotated[bool, typer.Option("--verbose", "-v", help="详细输出")] = False,
    quiet: Annotated[bool, typer.Option("--quiet", "-q", help="静默模式")] = False,
) -> None:
    """VCF处理器 - 独立脚本版本
    
    使用示例:
      python vcf_processor_standalone.py input.vcf output
      python vcf_processor_standalone.py input.vcf output --targets targets.txt
      python vcf_processor_standalone.py input.vcf output --targets targets.txt --batch-size 20000
      python vcf_processor_standalone.py input.vcf output --compress --verbose

    等价的模块化命令:
      python -m chip.vcf_processor.cli process input.vcf output
      python -m chip.vcf_processor.cli process input.vcf output --targets targets.txt
      python -m chip.vcf_processor.cli process input.vcf output --targets targets.txt --batch-size 20000

    输入文件格式:
      VCF文件: 标准VCF格式文件 (支持.vcf和.vcf.gz)
      目标文件: 每行一个变异ID，格式为 CHROM_POS_REF_ALT (可选)

    输出文件:
      output.genotype_codes.tsv: VCF格式基因型表格 (0/0, 0/1, 1/1, ./.)
      output.genotype_bases.tsv: 碱基序列基因型表格 (AA, AT, TT, NN)

    注意:
      此脚本直接使用 vcf_processor.py 中的VCFProcessor类，需要相应的依赖库。
      如果不指定目标文件，将处理VCF文件中的所有变异。
      输出文件使用TSV格式，便于后续数据分析和处理。
      缺失基因型在codes文件中显示为 ./. ，在bases文件中显示为 NN 。
    """
    try:
        # 设置日志
        from chip.vcf_processor.logging_config import setup_logging
        setup_logging(
            verbose=verbose and not quiet,
            quiet=quiet,
            log_file=None
        )
        
        # 使用VCFProcessor处理
        result = process_vcf(
            vcf_file=str(vcf_file),
            output_prefix=str(output_prefix),
            target_file=str(targets) if targets else None,
            miss_fmt=miss_fmt,
            gt_sep=gt_sep,
            batch_size=batch_size,
            compress_output=compress,
            verbose=verbose and not quiet,
            quiet=quiet
        )
        
        # 输出结果摘要
        if not quiet:
            typer.echo("\n" + "="*50)
            typer.echo("VCF处理完成")
            typer.echo("="*50)
            typer.echo(f"处理的变异数: {result['processed_variants']:,}")
            typer.echo(f"总变异数: {result['total_variants']:,}")
            typer.echo(f"处理时间: {result['elapsed_time']:.2f} 秒")
            typer.echo(f"处理速度: {result['processing_rate']:.1f} 变异/秒")
            
            if result['errors']:
                typer.echo(f"错误数: {len(result['errors'])}")
                for error in result['errors'][:3]:  # 显示前3个错误
                    typer.echo(f"  - {error}")
                if len(result['errors']) > 3:
                    typer.echo(f"  ... 还有 {len(result['errors']) - 3} 个错误")
            
            if result['success']:
                typer.echo("状态: ✓ 处理成功", color=typer.colors.GREEN)
            else:
                typer.echo("状态: ✗ 处理失败", color=typer.colors.RED)
            
            typer.echo("="*50)
        
        # 返回适当的退出码
        if not result['success']:
            raise typer.Exit(1)
        
    except Exception as e:
        typer.echo(f"错误: {e}", err=True, color=typer.colors.RED)
        raise typer.Exit(1)


if __name__ == '__main__':
    typer.run(main)