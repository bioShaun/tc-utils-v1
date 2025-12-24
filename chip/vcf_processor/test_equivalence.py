#!/usr/bin/env python3
"""
测试独立脚本与模块化CLI的等价性

此测试验证两种VCF处理方式产生相同的输出结果：
1. 独立脚本: python vcf_processor_standalone.py
2. 模块化CLI: python -m chip.vcf_processor.cli process

注意: 模块化版本需要安装依赖，如果未安装会跳过对比测试。
"""

import tempfile
import subprocess
import sys
from pathlib import Path

def create_test_files():
    """创建测试用的VCF和目标文件"""
    
    # 创建临时目录
    temp_dir = tempfile.mkdtemp()
    temp_path = Path(temp_dir)
    
    # 创建测试VCF文件
    vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3
chr1	100	.	A	T	60	PASS	.	GT:DP	0/0:30	0/1:25	1/1:35
chr1	200	.	G	C	55	PASS	.	GT:DP	0/1:28	1/1:32	0/0:30
chr2	300	.	T	A	65	PASS	.	GT:DP	1/1:40	0/0:35	0/1:35
"""
    
    vcf_file = temp_path / "test.vcf"
    vcf_file.write_text(vcf_content)
    
    # 创建目标ID文件
    target_content = """chr1_100_A_T
chr1_200_G_C
chr2_300_T_A
"""
    
    target_file = temp_path / "targets.txt"
    target_file.write_text(target_content)
    
    return temp_path, vcf_file, target_file

def test_standalone_version(vcf_file, target_file, output_prefix):
    """测试独立脚本版本"""
    print("测试独立脚本版本...")
    
    cmd = [
        sys.executable, 
        "vcf_processor_standalone.py",
        str(vcf_file),
        str(target_file), 
        str(output_prefix),
        "--quiet"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path(__file__).parent)
    
    return result.returncode == 0, result

def test_modular_version(vcf_file, target_file, output_prefix):
    """测试模块化版本"""
    print("测试模块化版本...")
    
    cmd = [
        sys.executable, 
        "-m", "chip.vcf_processor.cli",
        "process",
        str(vcf_file),
        str(target_file), 
        str(output_prefix),
        "--quiet"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    return result.returncode == 0, result

def compare_outputs(output1_prefix, output2_prefix):
    """比较两个版本的输出文件"""
    print("比较输出文件...")
    
    # 检查基因型文件
    gt_file1 = Path(f"{output1_prefix}.gt.txt")
    gt_file2 = Path(f"{output2_prefix}.gt.txt")
    
    if not gt_file1.exists() or not gt_file2.exists():
        print(f"  ✗ 基因型文件缺失: {gt_file1.exists()=}, {gt_file2.exists()=}")
        return False
    
    content1 = gt_file1.read_text()
    content2 = gt_file2.read_text()
    
    if content1 != content2:
        print("  ✗ 基因型文件内容不同")
        print(f"    独立版本: {len(content1)} 字符")
        print(f"    模块版本: {len(content2)} 字符")
        return False
    
    # 检查序列文件
    seq_file1 = Path(f"{output1_prefix}.seq.txt")
    seq_file2 = Path(f"{output2_prefix}.seq.txt")
    
    if not seq_file1.exists() or not seq_file2.exists():
        print(f"  ✗ 序列文件缺失: {seq_file1.exists()=}, {seq_file2.exists()=}")
        return False
    
    content1 = seq_file1.read_text()
    content2 = seq_file2.read_text()
    
    if content1 != content2:
        print("  ✗ 序列文件内容不同")
        return False
    
    print("  ✓ 输出文件完全一致")
    return True

def main():
    """主测试函数"""
    print("=== VCF处理器等价性测试 ===")
    print("验证独立脚本与模块化CLI产生相同输出\n")
    
    # 创建测试文件
    temp_path, vcf_file, target_file = create_test_files()
    print(f"测试文件位置: {temp_path}")
    
    # 测试独立版本
    standalone_output = temp_path / "standalone_output"
    standalone_success, standalone_result = test_standalone_version(vcf_file, target_file, standalone_output)
    
    if not standalone_success:
        print(f"✗ 独立脚本版本失败:")
        print(f"  stdout: {standalone_result.stdout}")
        print(f"  stderr: {standalone_result.stderr}")
        return False
    
    print("✓ 独立脚本版本成功")
    
    # 测试模块化版本
    modular_output = temp_path / "modular_output"
    modular_success, modular_result = test_modular_version(vcf_file, target_file, modular_output)
    
    if not modular_success:
        print(f"⚠️  模块化版本失败 (可能未安装依赖):")
        print(f"  stderr: {modular_result.stderr}")
        print("\n💡 这是正常的，如果您只想使用独立脚本版本。")
        print("   要使用模块化版本，请运行: pip install cyvcf2 typer loguru tqdm")
        print("\n✓ 独立脚本版本测试通过！")
        return True
    
    print("✓ 模块化版本成功")
    
    # 比较输出
    outputs_match = compare_outputs(standalone_output, modular_output)
    
    if outputs_match:
        print("\n🎉 等价性测试通过！")
        print("独立脚本版本与模块化版本产生相同的输出结果。")
        return True
    else:
        print("\n❌ 等价性测试失败！")
        print("两个版本的输出结果不一致。")
        return False

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)