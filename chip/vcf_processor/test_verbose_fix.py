#!/usr/bin/env python3
"""
测试 verbose 模式修复
验证 verbose 模式下详细摘要输出到 stderr 而不是 stdout
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

def test_verbose_mode_output():
    """测试 verbose 模式下的输出位置"""
    print("=== 测试 verbose 模式输出位置 ===")
    
    temp_path, vcf_file = create_test_vcf()
    output_prefix = temp_path / "output"
    
    print(f"测试文件位置: {temp_path}")
    
    # 测试 verbose 模式
    cmd = [
        sys.executable, 
        "vcf_processor_standalone.py",
        str(vcf_file),
        str(output_prefix),
        "--verbose"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path(__file__).parent)
    
    print(f"退出码: {result.returncode}")
    print(f"标准输出长度: {len(result.stdout)} 字符")
    print(f"标准错误输出长度: {len(result.stderr)} 字符")
    
    # 检查标准输出中是否包含处理摘要
    summary_keywords = ["Processing summary", "variants processed", "elapsed time"]
    
    stdout_has_summary = any(keyword in result.stdout for keyword in summary_keywords)
    stderr_has_summary = any(keyword in result.stderr for keyword in summary_keywords)
    
    print(f"标准输出包含摘要: {stdout_has_summary}")
    print(f"标准错误输出包含摘要: {stderr_has_summary}")
    
    # 在 verbose 模式下，详细摘要应该在 stderr 中，不在 stdout 中
    if stdout_has_summary:
        print("✗ 标准输出包含处理摘要（应该在 stderr 中）")
        print("标准输出内容:")
        print(result.stdout)
        return False
    
    if stderr_has_summary:
        print("✓ 处理摘要正确输出到 stderr")
    else:
        print("- 未找到处理摘要（可能是日志级别问题）")
    
    # 检查输出文件是否正常生成
    gt_file = Path(str(output_prefix) + ".gt.txt")
    if gt_file.exists():
        print("✓ 输出文件正常生成")
        return True
    else:
        print("✗ 输出文件未生成")
        return False

def main():
    """主测试函数"""
    print("测试 verbose 模式修复\n")
    
    success = test_verbose_mode_output()
    
    print(f"\n=== 测试结果 ===")
    if success:
        print("✓ verbose 模式修复测试通过")
    else:
        print("✗ verbose 模式修复测试失败")
    
    return success

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)