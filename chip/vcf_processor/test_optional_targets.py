#!/usr/bin/env python3
"""
测试可选目标文件功能
"""

import tempfile
import os
from pathlib import Path

def create_test_vcf():
    """创建测试VCF文件"""
    
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
chr2	400	.	C	G	50	PASS	.	GT:DP	0/0:25	0/1:30	1/1:30
chr3	500	.	A	G	70	PASS	.	GT:DP	0/1:35	1/1:45	0/0:40
"""
    
    vcf_file = temp_path / "test.vcf"
    vcf_file.write_text(vcf_content)
    
    # 创建目标ID文件（只包含部分变异）
    target_content = """chr1_100_A_T
chr2_300_T_A
chr3_500_A_G
"""
    
    target_file = temp_path / "targets.txt"
    target_file.write_text(target_content)
    
    return temp_path, vcf_file, target_file

def test_with_targets():
    """测试使用目标文件"""
    print("=== 测试使用目标文件 ===")
    
    temp_path, vcf_file, target_file = create_test_vcf()
    output_prefix = temp_path / "with_targets"
    
    # 导入独立处理器
    import sys
    sys.path.insert(0, str(Path(__file__).parent))
    
    from vcf_processor_standalone import SimpleVCFProcessor
    
    # 创建处理器（使用目标文件）
    processor = SimpleVCFProcessor(
        vcf_file=str(vcf_file),
        output_prefix=str(output_prefix),
        target_file=str(target_file),
        verbose=True
    )
    
    # 执行处理
    result = processor.process()
    
    print(f"处理结果:")
    print(f"  处理的变异数: {result['processed_variants']}")
    print(f"  总变异数: {result['total_variants']}")
    print(f"  成功: {result['success']}")
    
    # 检查输出文件
    gt_file = output_prefix.with_suffix('.gt.txt')
    if gt_file.exists():
        with open(gt_file, 'r') as f:
            lines = f.readlines()
            print(f"  输出行数: {len(lines)} (包含头部)")
            print("  前几行内容:")
            for line in lines[:4]:
                print(f"    {line.strip()}")
    
    return result['success'], result['processed_variants']

def test_without_targets():
    """测试不使用目标文件（处理所有变异）"""
    print("\n=== 测试不使用目标文件（处理所有变异）===")
    
    temp_path, vcf_file, target_file = create_test_vcf()
    output_prefix = temp_path / "without_targets"
    
    # 导入独立处理器
    import sys
    sys.path.insert(0, str(Path(__file__).parent))
    
    from vcf_processor_standalone import SimpleVCFProcessor
    
    # 创建处理器（不使用目标文件）
    processor = SimpleVCFProcessor(
        vcf_file=str(vcf_file),
        output_prefix=str(output_prefix),
        target_file=None,  # 不使用目标文件
        verbose=True
    )
    
    # 执行处理
    result = processor.process()
    
    print(f"处理结果:")
    print(f"  处理的变异数: {result['processed_variants']}")
    print(f"  总变异数: {result['total_variants']}")
    print(f"  成功: {result['success']}")
    
    # 检查输出文件
    gt_file = output_prefix.with_suffix('.gt.txt')
    if gt_file.exists():
        with open(gt_file, 'r') as f:
            lines = f.readlines()
            print(f"  输出行数: {len(lines)} (包含头部)")
            print("  前几行内容:")
            for line in lines[:6]:  # 显示更多行，因为处理了所有变异
                print(f"    {line.strip()}")
    
    return result['success'], result['processed_variants']

def main():
    """主测试函数"""
    print("测试可选目标文件功能\n")
    
    # 测试使用目标文件
    success1, count1 = test_with_targets()
    
    # 测试不使用目标文件
    success2, count2 = test_without_targets()
    
    # 比较结果
    print(f"\n=== 结果对比 ===")
    print(f"使用目标文件: 处理了 {count1} 个变异")
    print(f"不使用目标文件: 处理了 {count2} 个变异")
    
    if success1 and success2:
        if count2 > count1:
            print("✓ 测试通过！不使用目标文件时处理了更多变异（符合预期）")
            return True
        else:
            print("✗ 测试失败！不使用目标文件时应该处理更多变异")
            return False
    else:
        print("✗ 测试失败！处理过程中出现错误")
        return False

if __name__ == '__main__':
    success = main()
    exit(0 if success else 1)