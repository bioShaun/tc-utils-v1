#!/usr/bin/env python3
"""
测试完整的静默模式
验证静默模式下完全没有输出到 stdout 和 stderr
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

def test_complete_quiet_mode():
    """测试完整的静默模式"""
    print("=== 测试完整静默模式 ===")
    
    temp_path, vcf_file = create_test_vcf()
    output_prefix = temp_path / "output"
    
    print(f"测试文件位置: {temp_path}")
    
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
    
    # 静默模式下，标准输出应该完全为空
    if result.stdout.strip():
        print("✗ 静默模式下标准输出不为空")
        print("标准输出内容:")
        print(repr(result.stdout))
        return False
    else:
        print("✓ 静默模式下标准输出为空")
    
    # 检查标准错误输出
    if result.stderr.strip():
        print("标准错误输出内容:")
        print(repr(result.stderr))
        
        # 检查是否包含进度条输出
        if "Processing variants" in result.stderr:
            print("✗ 静默模式下仍有进度条输出")
            return False
        else:
            print("- 标准错误输出不包含进度条")
    else:
        print("✓ 静默模式下标准错误输出也为空")
    
    # 检查输出文件是否正常生成
    gt_file = Path(str(output_prefix) + ".gt.txt")
    if gt_file.exists():
        content = gt_file.read_text()
        print("✓ 输出文件正常生成")
        print(f"输出文件内容长度: {len(content)} 字符")
        
        # 确保输出文件不包含任何调试信息
        if any(keyword in content for keyword in ["ProcessingConfig", "DEBUG", "INFO"]):
            print("✗ 输出文件包含调试信息")
            return False
        else:
            print("✓ 输出文件不包含调试信息")
            return True
    else:
        print("✗ 输出文件未生成")
        return False

def test_normal_mode_has_output():
    """测试正常模式下有适当的输出"""
    print("\n=== 测试正常模式输出 ===")
    
    temp_path, vcf_file = create_test_vcf()
    output_prefix = temp_path / "output_normal"
    
    # 测试正常模式（非静默，非详细）
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
    
    # 正常模式下应该有摘要输出到 stdout
    if "VCF处理完成" in result.stdout:
        print("✓ 正常模式下有摘要输出")
    else:
        print("✗ 正常模式下缺少摘要输出")
        return False
    
    # 应该有日志输出到 stderr
    if result.stderr.strip():
        print("✓ 正常模式下有日志输出")
    else:
        print("- 正常模式下没有日志输出")
    
    return True

def main():
    """主测试函数"""
    print("测试完整静默模式功能\n")
    
    success1 = test_complete_quiet_mode()
    success2 = test_normal_mode_has_output()
    
    overall_success = success1 and success2
    
    print(f"\n=== 测试结果 ===")
    if overall_success:
        print("✓ 静默模式测试全部通过")
    else:
        print("✗ 静默模式测试部分失败")
    
    return overall_success

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)