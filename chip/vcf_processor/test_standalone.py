#!/usr/bin/env python3
"""
测试独立VCF处理脚本的功能
"""

import tempfile
import os
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
chr2	400	.	C	G	50	PASS	.	GT:DP	0/0:25	0/1:30	1/1:30
chr3	500	.	A	G	70	PASS	.	GT:DP	0/1:35	1/1:45	0/0:40
"""
    
    vcf_file = temp_path / "test.vcf"
    vcf_file.write_text(vcf_content)
    
    # 创建目标ID文件
    target_content = """chr1_100_A_T
chr1_200_G_C
chr2_300_T_A
chr3_500_A_G
"""
    
    target_file = temp_path / "targets.txt"
    target_file.write_text(target_content)
    
    return temp_path, vcf_file, target_file

def test_standalone_processor():
    """测试独立处理器"""
    print("创建测试文件...")
    temp_path, vcf_file, target_file = create_test_files()
    
    output_prefix = temp_path / "test_output"
    
    print(f"测试文件位置: {temp_path}")
    print(f"VCF文件: {vcf_file}")
    print(f"目标文件: {target_file}")
    print(f"输出前缀: {output_prefix}")
    
    # 导入并测试独立处理器
    import sys
    sys.path.insert(0, str(Path(__file__).parent))
    
    from vcf_processor_standalone import SimpleVCFProcessor
    
    # 创建处理器
    processor = SimpleVCFProcessor(
        vcf_file=str(vcf_file),
        target_file=str(target_file),
        output_prefix=str(output_prefix),
        miss_fmt="NN",
        gt_sep="",
        batch_size=1000,
        compress_output=False,
        verbose=True
    )
    
    # 执行处理
    print("\n开始处理...")
    result = processor.process()
    
    # 显示结果
    print("\n处理结果:")
    print(f"  处理的变异数: {result['processed_variants']}")
    print(f"  总变异数: {result['total_variants']}")
    print(f"  处理时间: {result['elapsed_time']:.2f} 秒")
    print(f"  处理速度: {result['processing_rate']:.1f} 变异/秒")
    print(f"  成功: {result['success']}")
    
    if result['errors']:
        print(f"  错误: {result['errors']}")
    
    # 检查输出文件
    gt_file = output_prefix.with_suffix('.gt.txt')
    seq_file = output_prefix.with_suffix('.seq.txt')
    
    print(f"\n输出文件:")
    if gt_file.exists():
        print(f"  ✓ 基因型文件: {gt_file}")
        print(f"    大小: {gt_file.stat().st_size} 字节")
        
        # 显示前几行内容
        with open(gt_file, 'r') as f:
            lines = f.readlines()[:5]
            print("    前几行内容:")
            for line in lines:
                print(f"      {line.strip()}")
    else:
        print(f"  ✗ 基因型文件未创建")
    
    if seq_file.exists():
        print(f"  ✓ 序列文件: {seq_file}")
        print(f"    大小: {seq_file.stat().st_size} 字节")
    else:
        print(f"  ✗ 序列文件未创建")
    
    print(f"\n测试完成！临时文件位置: {temp_path}")
    print("您可以查看生成的文件来验证结果。")
    
    return result['success']

if __name__ == '__main__':
    success = test_standalone_processor()
    if success:
        print("\n✓ 测试通过！")
    else:
        print("\n✗ 测试失败！")