#!/usr/bin/env python3
"""
测试独立脚本直接导入 vcf_processor.py
验证脚本能够成功导入并使用 VCFProcessor 类
"""

import tempfile
import subprocess
import sys
from pathlib import Path

def create_simple_vcf():
    """创建简单的测试VCF文件"""
    temp_dir = tempfile.mkdtemp()
    temp_path = Path(temp_dir)
    
    # 创建测试VCF文件
    vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2
chr1	100	.	A	T	60	PASS	.	GT	0/0	0/1
chr1	200	.	G	C	55	PASS	.	GT	0/1	1/1
"""
    
    vcf_file = temp_path / "test.vcf"
    vcf_file.write_text(vcf_content)
    
    return temp_path, vcf_file

def test_direct_import():
    """测试直接导入功能"""
    print("=== 测试直接导入 vcf_processor.py ===")
    
    temp_path, vcf_file = create_simple_vcf()
    output_prefix = temp_path / "output"
    
    print(f"测试文件位置: {temp_path}")
    
    # 测试处理
    cmd = [
        sys.executable, 
        "vcf_processor_standalone.py",
        str(vcf_file),
        str(output_prefix),
        "--verbose"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path(__file__).parent)
    
    print(f"退出码: {result.returncode}")
    print("标准错误输出:")
    print(result.stderr)
    
    # 检查是否使用了高性能版本
    if "Processing with VCFReader" in result.stderr and "processed 2 variants" in result.stderr:
        print("✓ 成功直接导入并使用了 vcf_processor.py")
        success = True
    elif "Batch processing completed: 2 total variants" in result.stderr:
        print("✓ 成功直接导入并使用了 vcf_processor.py")
        success = True
    elif "模块化版本不可用" in result.stderr:
        print("✗ 导入失败，回退到纯Python实现")
        success = False
    else:
        print("? 无法确定使用的版本")
        success = False
    
    # 检查输出文件
    gt_file = Path(str(output_prefix) + ".gt.txt")
    seq_file = Path(str(output_prefix) + ".seq.txt")
    if gt_file.exists() and seq_file.exists():
        print(f"✓ 输出文件已创建: {gt_file}, {seq_file}")
        # 检查文件内容
        if gt_file.stat().st_size > 0 and seq_file.stat().st_size > 0:
            print("✓ 输出文件包含数据")
        else:
            print("- 输出文件为空")
    else:
        print("- 输出文件未创建")
        success = False
    
    return success

def main():
    """主测试函数"""
    print("测试独立脚本直接导入功能\n")
    
    success = test_direct_import()
    
    print(f"\n=== 测试结果 ===")
    if success:
        print("✓ 直接导入测试通过")
    else:
        print("✗ 直接导入测试失败")
    
    return success

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)