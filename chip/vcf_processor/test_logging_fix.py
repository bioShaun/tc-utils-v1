#!/usr/bin/env python3
"""
测试日志输出修复
验证修复后的代码不会将日志输出混入数据输出
"""

import tempfile
import subprocess
import sys
from pathlib import Path

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
chr1	200	.	G	C	55	PASS	.	GT	0/1	1/1
"""
    
    vcf_file = temp_path / "test.vcf"
    vcf_file.write_text(vcf_content)
    
    return temp_path, vcf_file

def test_no_debug_output():
    """测试非 verbose 模式下不输出调试信息"""
    print("=== 测试非 verbose 模式下的日志输出 ===")
    
    temp_path, vcf_file = create_test_vcf()
    output_prefix = temp_path / "output"
    
    print(f"测试文件位置: {temp_path}")
    
    # 测试处理（非 verbose 模式）
    cmd = [
        sys.executable, 
        "vcf_processor_standalone.py",
        str(vcf_file),
        str(output_prefix)
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path(__file__).parent)
    
    print(f"退出码: {result.returncode}")
    print(f"标准输出长度: {len(result.stdout)} 字符")
    print(f"标准错误输出长度: {len(result.stderr)} 字符")
    
    # 检查标准错误输出中是否包含 ProcessingConfig
    if "ProcessingConfig" in result.stderr:
        print("✗ 标准错误输出包含 ProcessingConfig 调试信息")
        print("标准错误输出内容:")
        print(result.stderr)
        return False
    else:
        print("✓ 标准错误输出不包含 ProcessingConfig 调试信息")
    
    # 检查输出文件
    gt_file = Path(str(output_prefix) + ".gt.txt")
    
    if gt_file.exists():
        content = gt_file.read_text()
        print(f"✓ 输出文件已创建，内容长度: {len(content)} 字符")
        
        # 检查输出文件中是否包含 ProcessingConfig
        if "ProcessingConfig" in content:
            print("✗ 输出文件包含 ProcessingConfig 调试信息")
            print("输出文件内容:")
            print(content)
            return False
        else:
            print("✓ 输出文件不包含 ProcessingConfig 调试信息")
            print("输出文件内容:")
            print(content)
            return True
    else:
        print("✗ 输出文件未创建")
        return False

def test_verbose_mode():
    """测试 verbose 模式下的日志输出"""
    print("\n=== 测试 verbose 模式下的日志输出 ===")
    
    temp_path, vcf_file = create_test_vcf()
    output_prefix = temp_path / "output_verbose"
    
    # 测试处理（verbose 模式）
    cmd = [
        sys.executable, 
        "vcf_processor_standalone.py",
        str(vcf_file),
        str(output_prefix),
        "--verbose"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path(__file__).parent)
    
    print(f"退出码: {result.returncode}")
    print(f"标准错误输出长度: {len(result.stderr)} 字符")
    
    # 在 verbose 模式下，应该有调试输出到 stderr
    if "Configuration validated" in result.stderr:
        print("✓ verbose 模式下正确输出调试信息到 stderr")
    else:
        print("- verbose 模式下没有找到预期的调试信息")
    
    # 检查输出文件
    gt_file = Path(str(output_prefix) + ".gt.txt")
    
    if gt_file.exists():
        content = gt_file.read_text()
        
        # 即使在 verbose 模式下，输出文件也不应该包含调试信息
        if "ProcessingConfig" in content:
            print("✗ 输出文件包含 ProcessingConfig 调试信息")
            return False
        else:
            print("✓ 输出文件不包含 ProcessingConfig 调试信息")
            return True
    else:
        print("✗ 输出文件未创建")
        return False

def main():
    """主测试函数"""
    print("测试日志输出修复\n")
    
    success1 = test_no_debug_output()
    success2 = test_verbose_mode()
    
    overall_success = success1 and success2
    
    print(f"\n=== 测试结果 ===")
    if overall_success:
        print("✓ 日志输出修复测试全部通过")
    else:
        print("✗ 日志输出修复测试部分失败")
    
    return overall_success

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)