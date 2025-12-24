#!/usr/bin/env python3
"""
测试输出格式
验证.gt.txt输出原始VCF格式，.seq.txt输出序列格式
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
chr1	300	.	C	.	50	PASS	.	GT	0/0	./.
"""
    
    vcf_file = temp_path / "test.vcf"
    vcf_file.write_text(vcf_content)
    
    return temp_path, vcf_file

def test_output_formats():
    """测试输出格式"""
    print("=== 测试输出格式 ===")
    
    temp_path, vcf_file = create_test_vcf()
    output_prefix = temp_path / "output"
    
    print(f"测试文件位置: {temp_path}")
    
    # 测试处理
    cmd = [
        sys.executable, 
        "vcf_processor_standalone.py",
        str(vcf_file),
        str(output_prefix),
        "--quiet"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path(__file__).parent)
    
    if result.returncode != 0:
        print(f"✗ 处理失败，退出码: {result.returncode}")
        print(f"错误输出: {result.stderr}")
        return False
    
    # 检查输出文件
    gt_file = Path(str(output_prefix) + ".gt.txt")
    seq_file = Path(str(output_prefix) + ".seq.txt")
    
    if not gt_file.exists():
        print("✗ .gt.txt 文件未生成")
        return False
    
    if not seq_file.exists():
        print("✗ .seq.txt 文件未生成")
        return False
    
    # 检查 .gt.txt 文件格式（应该是原始VCF格式）
    gt_content = gt_file.read_text()
    print("GT文件内容:")
    print(gt_content)
    
    # 检查是否包含VCF格式的基因型
    if "0/0" in gt_content and "0/1" in gt_content and "1/1" in gt_content:
        print("✓ .gt.txt 文件包含原始VCF基因型格式")
    else:
        print("✗ .gt.txt 文件不包含预期的VCF基因型格式")
        return False
    
    # 检查 .seq.txt 文件格式（应该是序列格式）
    seq_content = seq_file.read_text()
    print("\nSEQ文件内容:")
    print(seq_content)
    
    # 检查是否包含序列格式的基因型
    if "AA" in seq_content and ("AT" in seq_content or "GC" in seq_content) and "CC" in seq_content:
        print("✓ .seq.txt 文件包含序列基因型格式")
    else:
        print("✗ .seq.txt 文件不包含预期的序列基因型格式")
        print(f"实际内容: {seq_content}")
        return False
    
    # 检查缺失基因型处理
    if "NN" in gt_content and "NN" in seq_content:
        print("✓ 缺失基因型正确处理为 NN")
    else:
        print("✗ 缺失基因型处理不正确")
        return False
    
    return True

def main():
    """主测试函数"""
    print("测试VCF处理器输出格式\n")
    
    success = test_output_formats()
    
    print(f"\n=== 测试结果 ===")
    if success:
        print("✓ 输出格式测试通过")
        print("  - .gt.txt: 原始VCF基因型格式 (0/0, 0/1, 1/1)")
        print("  - .seq.txt: 序列基因型格式 (AA, AT, TT)")
        print("  - 缺失基因型: NN")
    else:
        print("✗ 输出格式测试失败")
    
    return success

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)