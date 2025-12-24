#!/usr/bin/env python3
"""
测试标准输出污染问题
检查是否有任何调试信息泄漏到标准输出
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
chr1	200	.	G		55	PASS	.	GT	0/1	1/1
chr1	300	.	C	A,G	50	PASS	.	GT	0/1	1/2
"""
    
    vcf_file = temp_path / "test.vcf"
    vcf_file.write_text(vcf_content)
    
    return temp_path, vcf_file

def test_stdout_only():
    """测试只捕获标准输出，看是否有调试信息"""
    print("=== 测试标准输出污染 ===")
    
    temp_path, vcf_file = create_test_vcf()
    output_prefix = temp_path / "output"
    
    print(f"测试文件位置: {temp_path}")
    
    # 测试处理，只捕获标准输出
    cmd = [
        sys.executable, 
        "vcf_processor_standalone.py",
        str(vcf_file),
        str(output_prefix)
    ]
    
    # 分别捕获 stdout 和 stderr
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path(__file__).parent)
    
    print(f"退出码: {result.returncode}")
    print(f"标准输出内容:")
    print("=" * 40)
    print(repr(result.stdout))
    print("=" * 40)
    print(f"标准错误输出内容:")
    print("=" * 40)
    print(repr(result.stderr))
    print("=" * 40)
    
    # 检查标准输出中是否包含任何调试信息
    debug_keywords = ["ProcessingConfig", "Configuration validated", "DEBUG", "loguru"]
    
    stdout_contaminated = False
    for keyword in debug_keywords:
        if keyword in result.stdout:
            print(f"✗ 标准输出包含调试关键词: {keyword}")
            stdout_contaminated = True
    
    if not stdout_contaminated:
        print("✓ 标准输出不包含调试信息")
    
    return not stdout_contaminated

def test_quiet_mode():
    """测试静默模式"""
    print("\n=== 测试静默模式 ===")
    
    temp_path, vcf_file = create_test_vcf()
    output_prefix = temp_path / "output_quiet"
    
    # 测试静默模式
    cmd = [
        sys.executable, 
        "vcf_processor_standalone.py",
        str(vcf_file),
        str(output_prefix),
        "--quiet"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path(__file__).parent)
    
    print(f"退出码: {result.returncode}")
    print(f"标准输出长度: {len(result.stdout)} 字符")
    print(f"标准错误输出长度: {len(result.stderr)} 字符")
    
    # 静默模式下，标准输出应该为空
    if result.stdout.strip():
        print("✗ 静默模式下标准输出不为空")
        print("标准输出内容:")
        print(repr(result.stdout))
        return False
    else:
        print("✓ 静默模式下标准输出为空")
    
    # 检查是否有警告级别的输出到 stderr
    if result.stderr.strip():
        print(f"标准错误输出: {repr(result.stderr)}")
    
    return True

def test_with_empty_alt():
    """测试包含空ALT字段的VCF文件"""
    print("\n=== 测试空ALT字段处理 ===")
    
    temp_path, vcf_file = create_test_vcf()
    output_prefix = temp_path / "output_empty_alt"
    
    # 测试处理
    cmd = [
        sys.executable, 
        "vcf_processor_standalone.py",
        str(vcf_file),
        str(output_prefix)
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path(__file__).parent)
    
    print(f"退出码: {result.returncode}")
    
    # 检查输出文件
    gt_file = Path(str(output_prefix) + ".gt.txt")
    
    if gt_file.exists():
        content = gt_file.read_text()
        print("输出文件内容:")
        print(content)
        
        # 检查空ALT字段是否正确显示为 "."
        lines = content.strip().split('\n')
        if len(lines) > 2:  # 跳过标题行
            second_line = lines[2]  # chr1 200行应该有空ALT
            if '\t.\t' in second_line:
                print("✓ 空ALT字段正确显示为 '.'")
                return True
            else:
                print("✗ 空ALT字段显示不正确")
                return False
    
    return False

def main():
    """主测试函数"""
    print("测试标准输出污染问题\n")
    
    success1 = test_stdout_only()
    success2 = test_quiet_mode()
    success3 = test_with_empty_alt()
    
    overall_success = success1 and success2 and success3
    
    print(f"\n=== 测试结果 ===")
    if overall_success:
        print("✓ 所有测试通过，无标准输出污染")
    else:
        print("✗ 发现标准输出污染问题")
    
    return overall_success

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)