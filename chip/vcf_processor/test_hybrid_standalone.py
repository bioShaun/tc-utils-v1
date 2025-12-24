#!/usr/bin/env python3
"""
测试混合独立脚本功能
验证脚本能够正确选择高性能模块或回退到纯Python实现
"""

import tempfile
import subprocess
import sys
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

def test_standalone_script():
    """测试独立脚本的命令行接口"""
    print("=== 测试混合独立脚本 ===")
    
    temp_path, vcf_file, target_file = create_test_vcf()
    output_prefix = temp_path / "test_output"
    
    print(f"测试文件位置: {temp_path}")
    
    # 测试处理所有变异
    print("\n--- 测试处理所有变异 ---")
    cmd_all = [
        sys.executable, 
        "vcf_processor_standalone.py",
        str(vcf_file),
        str(output_prefix) + "_all",
        "--verbose"
    ]
    
    result_all = subprocess.run(cmd_all, capture_output=True, text=True, cwd=Path(__file__).parent)
    
    print(f"退出码: {result_all.returncode}")
    if result_all.stdout:
        print("标准输出:")
        print(result_all.stdout)
    if result_all.stderr:
        print("标准错误:")
        print(result_all.stderr)
    
    # 测试处理目标变异
    print("\n--- 测试处理目标变异 ---")
    cmd_targets = [
        sys.executable, 
        "vcf_processor_standalone.py",
        str(vcf_file),
        str(output_prefix) + "_targets",
        "--targets", str(target_file),
        "--verbose"
    ]
    
    result_targets = subprocess.run(cmd_targets, capture_output=True, text=True, cwd=Path(__file__).parent)
    
    print(f"退出码: {result_targets.returncode}")
    if result_targets.stdout:
        print("标准输出:")
        print(result_targets.stdout)
    if result_targets.stderr:
        print("标准错误:")
        print(result_targets.stderr)
    
    # 检查输出文件
    print("\n--- 检查输出文件 ---")
    
    # 检查处理所有变异的输出
    gt_file_all = Path(str(output_prefix) + "_all.gt.txt")
    if gt_file_all.exists():
        with open(gt_file_all, 'r') as f:
            lines_all = f.readlines()
            print(f"处理所有变异: {len(lines_all)-1} 个变异 (不含头部)")
    else:
        print("✗ 处理所有变异的输出文件未找到")
    
    # 检查处理目标变异的输出
    gt_file_targets = Path(str(output_prefix) + "_targets.gt.txt")
    if gt_file_targets.exists():
        with open(gt_file_targets, 'r') as f:
            lines_targets = f.readlines()
            print(f"处理目标变异: {len(lines_targets)-1} 个变异 (不含头部)")
    else:
        print("✗ 处理目标变异的输出文件未找到")
    
    # 验证结果
    success = (result_all.returncode == 0 and result_targets.returncode == 0 and
               gt_file_all.exists() and gt_file_targets.exists())
    
    if success:
        print("\n✓ 混合独立脚本测试通过！")
        
        # 显示使用的处理器类型
        if "使用高性能模块化版本" in result_all.stderr:
            print("  - 成功使用了高性能模块化版本")
        elif "使用纯Python实现" in result_all.stderr:
            print("  - 回退到纯Python实现 (正常，如果依赖未安装)")
        else:
            print("  - 处理器类型检测: 未明确显示")
            
    else:
        print("\n✗ 混合独立脚本测试失败！")
    
    return success

def test_help_message():
    """测试帮助信息"""
    print("\n=== 测试帮助信息 ===")
    
    cmd = [sys.executable, "vcf_processor_standalone.py", "--help"]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path(__file__).parent)
    
    print("帮助信息:")
    print(result.stdout)
    
    # 检查关键信息是否存在
    help_text = result.stdout
    checks = [
        "智能独立脚本版本" in help_text,
        "智能处理模式" in help_text,
        "--targets" in help_text,
        "可选" in help_text
    ]
    
    if all(checks):
        print("✓ 帮助信息包含所有必要内容")
        return True
    else:
        print("✗ 帮助信息缺少某些内容")
        return False

def main():
    """主测试函数"""
    print("测试混合独立脚本功能\n")
    
    # 测试帮助信息
    help_success = test_help_message()
    
    # 测试脚本功能
    script_success = test_standalone_script()
    
    # 总结
    print(f"\n=== 测试总结 ===")
    print(f"帮助信息测试: {'✓ 通过' if help_success else '✗ 失败'}")
    print(f"脚本功能测试: {'✓ 通过' if script_success else '✗ 失败'}")
    
    overall_success = help_success and script_success
    print(f"总体结果: {'✓ 全部通过' if overall_success else '✗ 有测试失败'}")
    
    return overall_success

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)