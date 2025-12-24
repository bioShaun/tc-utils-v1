#!/usr/bin/env python3
"""
调试静默模式问题
"""

import tempfile
import sys
from pathlib import Path

# 添加项目根目录到Python路径以支持绝对导入
current_dir = Path(__file__).parent
project_root = current_dir.parent.parent  # 从 chip/vcf_processor 到项目根目录
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from chip.vcf_processor.config import ProcessingConfig
from chip.vcf_processor.vcf_processor import VCFProcessor

def create_test_vcf():
    """创建测试VCF文件"""
    temp_dir = tempfile.mkdtemp()
    temp_path = Path(temp_dir)
    
    # 创建测试VCF文件
    vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2
chr1	100	.	A	T	60	PASS	.	GT	0/0	0/1
chr1	200	.	G		55	PASS	.	GT	0/1	1/1
chr1	300	.	C	A,G	50	PASS	.	GT	0/1	1/2
"""
    
    vcf_file = temp_path / "test.vcf"
    vcf_file.write_text(vcf_content)
    
    return temp_path, vcf_file

def test_direct_quiet_mode():
    """直接测试静默模式"""
    print("=== 直接测试静默模式 ===")
    
    temp_path, vcf_file = create_test_vcf()
    output_prefix = temp_path / "output"
    
    print(f"测试文件位置: {temp_path}")
    
    # 创建配置
    config = ProcessingConfig(
        vcf_file=vcf_file,
        output_file=output_prefix,
        target_id_file=None,
        miss_fmt="NN",
        gt_sep="",
        batch_size=10000,
        compress_output=False,
        verbose=False,
        quiet=True,  # 静默模式
        threads=1
    )
    
    print(f"配置 quiet: {config.quiet}")
    print(f"配置 verbose: {config.verbose}")
    
    # 执行处理
    processor = VCFProcessor(config)
    result = processor.process()
    
    print(f"处理结果: {result.processed_variants} 个变异")
    
    return True

if __name__ == '__main__':
    test_direct_quiet_mode()