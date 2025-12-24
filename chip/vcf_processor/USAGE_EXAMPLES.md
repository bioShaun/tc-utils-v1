# VCF处理器使用示例

本文档展示VCF处理器两种使用方式的具体示例。

## 示例数据准备

首先创建测试数据：

```bash
# 创建示例VCF文件
cat > example.vcf << 'EOF'
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3
chr1	100	.	A	T	60	PASS	.	GT	0/0	0/1	1/1
chr1	200	.	G	C	55	PASS	.	GT	0/1	1/1	0/0
chr2	300	.	T	A	65	PASS	.	GT	1/1	0/0	0/1
chr2	400	.	C	G	50	PASS	.	GT	0/0	0/1	1/1
EOF

# 创建目标变异ID文件（可选）
cat > targets.txt << 'EOF'
chr1_100_A_T
chr1_200_G_C
chr2_300_T_A
EOF
```

## 方式1: 独立脚本 (推荐快速开始)

### 基本使用

```bash
# 处理所有变异（最简单的使用方式）
python vcf_processor_standalone.py example.vcf results

# 处理特定目标变异
python vcf_processor_standalone.py example.vcf results --targets targets.txt

# 输出:
# ==================================================
# VCF处理完成
# ==================================================
# 处理的变异数: 3 (使用目标文件) 或 4 (处理所有变异)
# 总变异数: 4
# 处理时间: 0.01 秒
# 处理速度: 300.0 变异/秒
# 状态: ✓ 处理成功
# ==================================================
```

### 自定义参数

```bash
# 使用自定义参数处理所有变异
python vcf_processor_standalone.py example.vcf results \
    --miss-fmt "./." \
    --gt-sep "/" \
    --batch-size 1000 \
    --compress \
    --verbose

# 使用自定义参数处理目标变异
python vcf_processor_standalone.py example.vcf results \
    --targets targets.txt \
    --miss-fmt "./." \
    --gt-sep "/" \
    --batch-size 1000 \
    --compress \
    --verbose
```

### 静默模式

```bash
# 静默处理所有变异，适合脚本自动化
python vcf_processor_standalone.py example.vcf results --quiet
echo "Exit code: $?"

# 静默处理目标变异
python vcf_processor_standalone.py example.vcf results --targets targets.txt --quiet
echo "Exit code: $?"
```

## 方式2: 模块化命令 (推荐生产环境)

### 前提条件

```bash
# 安装依赖 (仅需一次)
pip install cyvcf2 typer loguru tqdm pandas
```

### 基本使用

```bash
# 与独立脚本完全相同的命令格式
# 处理所有变异
python -m chip.vcf_processor.cli process example.vcf results

# 处理目标变异
python -m chip.vcf_processor.cli process example.vcf results --targets targets.txt

# 输出格式与独立脚本相同，但处理速度更快
```

### 高级功能

```bash
# 使用多线程处理所有变异 (独立脚本不支持)
python -m chip.vcf_processor.cli process example.vcf results \
    --threads 4 \
    --batch-size 20000 \
    --compress \
    --verbose

# 使用多线程处理目标变异
python -m chip.vcf_processor.cli process example.vcf results \
    --targets targets.txt \
    --threads 4 \
    --batch-size 20000 \
    --compress \
    --verbose

# 使用配置文件 (独立脚本不支持)
python -m chip.vcf_processor.cli create-config config.json
python -m chip.vcf_processor.cli process example.vcf results --config config.json
```

## 输出文件对比

两种方式产生完全相同的输出文件：

### 基因型文件 (results.gt.txt)

```
CHROM	POS	REF	ALT	sample1	sample2	sample3
chr1	100	A	T	AA	AT	TT
chr1	200	G	C	GC	CC	GG
chr2	300	T	A	AA	TT	TA
```

### 序列文件 (results.seq.txt)

```
CHROM	POS	REF	ALT	sample1	sample2	sample3
chr1	100	A	T	AA	AT	TT
chr1	200	G	C	GC	CC	GG
chr2	300	T	A	AA	TT	TA
```

## 性能对比示例

### 小文件处理 (<1MB)

```bash
# 独立脚本 - 足够快
time python vcf_processor_standalone.py small.vcf targets.txt output1
# 输出: real 0m0.123s

# 模块化版本 - 略快
time python -m chip.vcf_processor.cli process small.vcf targets.txt output2
# 输出: real 0m0.089s
```

### 大文件处理 (>100MB)

```bash
# 独立脚本 - 较慢但可用
time python vcf_processor_standalone.py large.vcf.gz targets.txt output1 \
    --batch-size 50000 --compress --quiet
# 输出: real 2m15.456s

# 模块化版本 - 显著更快
time python -m chip.vcf_processor.cli process large.vcf.gz targets.txt output2 \
    --batch-size 50000 --threads 4 --compress --quiet
# 输出: real 0m45.123s
```

## 错误处理示例

### 文件不存在

```bash
# 独立脚本
python vcf_processor_standalone.py missing.vcf targets.txt output
# 输出: 错误: VCF文件不存在: missing.vcf

# 模块化版本
python -m chip.vcf_processor.cli process missing.vcf targets.txt output
# 输出: 类似的错误信息，但格式可能略有不同
```

### 格式错误

```bash
# 两种方式都会检测并报告VCF格式错误
python vcf_processor_standalone.py invalid.vcf targets.txt output
python -m chip.vcf_processor.cli process invalid.vcf targets.txt output
```

## 脚本集成示例

### Bash脚本集成

```bash
#!/bin/bash
# process_vcf.sh

VCF_FILE="$1"
TARGETS="$2"
OUTPUT="$3"

# 检查是否安装了模块化版本
if python -c "import chip.vcf_processor.cli" 2>/dev/null; then
    echo "使用模块化版本 (高性能)"
    python -m chip.vcf_processor.cli process "$VCF_FILE" "$TARGETS" "$OUTPUT" --quiet
else
    echo "使用独立脚本版本 (零依赖)"
    python vcf_processor_standalone.py "$VCF_FILE" "$TARGETS" "$OUTPUT" --quiet
fi

echo "处理完成，退出码: $?"
```

### Python脚本集成

```python
#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path

def process_vcf(vcf_file, targets_file, output_prefix):
    """使用可用的VCF处理器版本"""
    
    # 尝试使用模块化版本
    try:
        import chip.vcf_processor.cli
        cmd = [
            sys.executable, "-m", "chip.vcf_processor.cli", "process",
            str(vcf_file), str(targets_file), str(output_prefix), "--quiet"
        ]
        print("使用模块化版本 (高性能)")
    except ImportError:
        # 回退到独立脚本
        script_path = Path(__file__).parent / "vcf_processor_standalone.py"
        cmd = [
            sys.executable, str(script_path),
            str(vcf_file), str(targets_file), str(output_prefix), "--quiet"
        ]
        print("使用独立脚本版本 (零依赖)")
    
    # 执行处理
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode == 0:
        print("✓ 处理成功")
        return True
    else:
        print(f"✗ 处理失败: {result.stderr}")
        return False

if __name__ == "__main__":
    success = process_vcf("example.vcf", "targets.txt", "results")
    sys.exit(0 if success else 1)
```

## 选择建议

### 使用独立脚本的场景:
- 🎯 快速测试和原型开发
- 🎯 无法安装额外依赖的环境
- 🎯 一次性数据处理任务
- 🎯 教学和演示
- 🎯 小到中等规模的文件 (<100MB)

### 使用模块化版本的场景:
- 🚀 生产环境和自动化流程
- 🚀 大规模数据处理 (>100MB)
- 🚀 需要高性能和多线程
- 🚀 需要配置文件管理
- 🚀 长期维护的项目

### 混合使用策略:
1. **开发阶段**: 使用独立脚本快速验证功能
2. **测试阶段**: 使用两种方式验证结果一致性
3. **生产阶段**: 使用模块化版本获得最佳性能

两种方式的命令参数完全兼容，可以随时无缝切换！