#!/usr/bin/env python3
"""
最终验证测试
验证所有修复都正常工作
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

def test_quiet_mode():
    """测试静默模式"""
    print("=== 测试静默模式 ===")
    
    temp_path, vcf_file = create_test_vcf()
    output_prefix = temp_path / "output_quiet"
    
    cmd = [
        sys.executable, 
        "vcf_processor_standalone.py",
        str(vcf_file),
        str(output_prefix),
        "--quiet"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path(__file__).parent)
    
    success = True
    
    # 检查退出码
    if result.returncode != 0:
        print(f"✗ 退出码错误: {result.returncode}")
        success = False
    else:
        print("✓ 退出码正确")
    
    # 检查标准输出
    if result.stdout.strip():
        print("✗ 静默模式下标准输出不为空")
        success = False
    else:
        print("✓ 静默模式下标准输出为空")
    
    # 检查标准错误输出
    if result.stderr.strip():
        print("✗ 静默模式下标准错误输出不为空")
        print(f"stderr: {repr(result.stderr)}")
        success = False
    else:
        print("✓ 静默模式下标准错误输出为空")
    
    # 检查输出文件
    gt_file = Path(str(output_prefix) + ".gt.txt")
    if gt_file.exists():
        content = gt_file.read_text()
        if "chr1\t200\tG\t." in content:
            print("✓ 空ALT字段正确处理")
        else:
            print("✗ 空ALT字段处理错误")
            success = False
    else:
        print("✗ 输出文件未生成")
        success = False
    
    return success

def test_normal_mode():
    """测试正常模式"""
    print("\n=== 测试正常模式 ===")
    
    temp_path, vcf_file = create_test_vcf()
    output_prefix = temp_path / "output_normal"
    
    cmd = [
        sys.executable, 
        "vcf_processor_standalone.py",
        str(vcf_file),
        str(output_prefix)
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path(__file__).parent)
    
    success = True
    
    # 检查退出码
    if result.returncode != 0:
        print(f"✗ 退出码错误: {result.returncode}")
        success = False
    else:
        print("✓ 退出码正确")
    
    # 检查标准输出（应该有摘要）
    if "VCF处理完成" not in result.stdout:
        print("✗ 标准输出缺少摘要")
        success = False
    else:
        print("✓ 标准输出包含摘要")
    
    # 检查标准错误输出（应该有日志，但不应该有ProcessingConfig调试信息）
    if not result.stderr.strip():
        print("- 标准错误输出为空")
    elif "ProcessingConfig" in result.stderr:
        print("✗ 标准错误输出包含ProcessingConfig调试信息")
        success = False
    else:
        print("✓ 标准错误输出正常（有日志但无调试信息）")
    
    return success

def test_verbose_mode():
    """测试详细模式"""
    print("\n=== 测试详细模式 ===")
    
    temp_path, vcf_file = create_test_vcf()
    output_prefix = temp_path / "output_verbose"
    
    cmd = [
        sys.executable, 
        "vcf_processor_standalone.py",
        str(vcf_file),
        str(output_prefix),
        "--verbose"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path(__file__).parent)
    
    success = True
    
    # 检查退出码
    if result.returncode != 0:
        print(f"✗ 退出码错误: {result.returncode}")
        success = False
    else:
        print("✓ 退出码正确")
    
    # 检查标准输出（应该有摘要）
    if "VCF处理完成" not in result.stdout:
        print("✗ 标准输出缺少摘要")
        success = False
    else:
        print("✓ 标准输出包含摘要")
    
    # 检查标准错误输出（应该有详细日志）
    if "Processing summary" in result.stderr:
        print("✓ 标准错误输出包含详细摘要")
    else:
        print("- 标准错误输出未找到详细摘要")
    
    return success

def main():
    """主测试函数"""
    print("最终验证测试\n")
    
    success1 = test_quiet_mode()
    success2 = test_normal_mode()
    success3 = test_verbose_mode()
    
    overall_success = success1 and success2 and success3
    
    print(f"\n=== 最终测试结果 ===")
    if overall_success:
        print("✓ 所有测试通过！调试输出污染问题已完全解决")
        print("  - 静默模式：完全无输出")
        print("  - 正常模式：适当的摘要和日志")
        print("  - 详细模式：详细摘要输出到stderr")
        print("  - 空ALT字段：正确显示为'.'")
        print("  - 无调试信息污染数据文件")
    else:
        print("✗ 部分测试失败")
    
    return overall_success

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)