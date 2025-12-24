# VCF处理器 - 独立脚本使用说明

## 概述

`vcf_processor_standalone.py` 是 `python -m chip.vcf_processor.cli process` 命令的**智能独立版本**。它会自动检测并使用最佳可用实现：

1. **优先使用高性能模块化版本** (如果cyvcf2等依赖可用)  
2. **自动回退到纯Python实现** (如果依赖不可用)

这种设计让您无需担心环境配置，脚本会自动适应您的系统！

## 智能适应特性

### 🚀 高性能模式 (自动检测)
- 当检测到 `cyvcf2`, `typer`, `loguru` 等依赖时自动启用
- 使用与模块化版本相同的高性能C扩展
- 处理速度可提升数倍

### 🛡️ 兼容模式 (自动回退)  
- 当依赖不可用时自动启用
- 仅使用Python标准库，零依赖
- 确保在任何环境下都能运行

## 使用方式

### 统一命令格式 (自动适应)
```bash
# 处理所有变异
python vcf_processor_standalone.py input.vcf output

# 处理特定目标变异  
python vcf_processor_standalone.py input.vcf output --targets targets.txt
```

脚本会自动选择最佳实现：
- ✅ **有依赖**: 自动使用高性能模块化版本
- ✅ **无依赖**: 自动回退到纯Python实现
- ✅ **任何环境**: 都能正常工作

### 与模块化命令的对比
```bash
# 独立脚本 (智能适应)
python vcf_processor_standalone.py input.vcf output --targets targets.txt

# 模块化命令 (需要安装)
python -m chip.vcf_processor.cli process input.vcf output --targets targets.txt
```

## 性能对比

| 环境 | 使用的实现 | 性能 | 依赖要求 | 自动检测 |
|------|------------|------|----------|----------|
| **完整环境** | 高性能模块化版本 | ⭐⭐⭐⭐⭐ | cyvcf2, typer等 | ✅ 自动启用 |
| **受限环境** | 纯Python实现 | ⭐⭐⭐ | 仅Python标准库 | ✅ 自动回退 |
| **任何环境** | 智能选择最佳 | ⭐⭐⭐→⭐⭐⭐⭐⭐ | 自适应 | ✅ 无需配置 |

## 智能脚本特点

- 🧠 **智能适应**: 自动检测环境并选择最佳实现
- 🚀 **性能优化**: 优先使用高性能C扩展 (如果可用)
- 🛡️ **兼容保证**: 在任何环境下都能运行
- ✅ **即用即走**: 下载脚本即可直接运行，无需配置
- ✅ **功能完整**: 支持VCF文件处理、基因型转换、批处理等核心功能
- ✅ **格式兼容**: 支持压缩和非压缩的VCF文件
- ✅ **中文友好**: 完整的中文提示和错误信息

## 快速开始

### 1. 下载脚本

```bash
# 下载脚本文件
wget vcf_processor_standalone.py
# 或者直接复制脚本内容到本地文件
```

### 2. 基本使用

```bash
# 最简单的使用方式
python vcf_processor_standalone.py input.vcf targets.txt output

# 查看帮助信息
python vcf_processor_standalone.py --help
```

## 详细使用说明

### 命令格式

```bash
python vcf_processor_standalone.py <VCF文件> <目标文件> <输出前缀> [选项]
```

### 必需参数

1. **VCF文件**: 输入的VCF格式文件路径
   - 支持 `.vcf` 和 `.vcf.gz` 格式
   - 必须包含标准VCF头部信息

2. **目标文件**: 包含目标变异ID的文本文件
   - 每行一个变异ID
   - 格式: `CHROM_POS_REF_ALT`
   - 例如: `chr1_100_A_T`

3. **输出前缀**: 输出文件的前缀名称
   - 会生成 `前缀.gt.txt` (基因型表) 和 `前缀.seq.txt` (序列表)

### 可选参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--miss-fmt` | `NN` | 缺失基因型的表示格式 |
| `--gt-sep` | `""` | 基因型等位基因间的分隔符 |
| `--batch-size` | `10000` | 批处理大小，影响内存使用和处理速度 |
| `--compress` | `False` | 是否压缩输出文件 (.gz格式) |
| `--verbose, -v` | `False` | 显示详细处理信息 |
| `--quiet, -q` | `False` | 静默模式，只显示错误信息 |

## 使用示例

### 示例1: 基本处理 (两种方式等价)

**独立脚本方式:**
```bash
python vcf_processor_standalone.py genotypes.vcf targets.txt results
```

**模块化方式:**
```bash
python -m chip.vcf_processor.cli process genotypes.vcf targets.txt results
```

**输入文件:**
- `genotypes.vcf`: VCF格式的基因型文件
- `targets.txt`: 目标变异ID列表

**输出文件:**
- `results.gt.txt`: 基因型表格
- `results.seq.txt`: 序列表格

### 示例2: 自定义参数 (两种方式等价)

**独立脚本方式:**
```bash
python vcf_processor_standalone.py genotypes.vcf.gz targets.txt results \
    --miss-fmt "./." \
    --gt-sep "/" \
    --batch-size 20000 \
    --compress \
    --verbose
```

**模块化方式:**
```bash
python -m chip.vcf_processor.cli process genotypes.vcf.gz targets.txt results \
    --miss-fmt "./." \
    --gt-sep "/" \
    --batch-size 20000 \
    --compress \
    --verbose
```

### 示例3: 大文件处理 (推荐使用模块化版本)

**独立脚本方式:**
```bash
# 处理大文件，使用大批处理大小和压缩输出
python vcf_processor_standalone.py large_file.vcf.gz targets.txt large_results \
    --batch-size 50000 \
    --compress \
    --verbose
```

**模块化方式 (推荐):**
```bash
# 模块化版本处理大文件更高效
python -m chip.vcf_processor.cli process large_file.vcf.gz targets.txt large_results \
    --batch-size 50000 \
    --threads 8 \
    --compress \
    --verbose
```

### 示例4: 静默处理 (两种方式等价)

**独立脚本方式:**
```bash
# 静默模式，适合脚本自动化
python vcf_processor_standalone.py input.vcf targets.txt output --quiet
```

**模块化方式:**
```bash
# 静默模式，适合脚本自动化
python -m chip.vcf_processor.cli process input.vcf targets.txt output --quiet
```

## 输入文件格式

### VCF文件格式

标准VCF格式文件，必须包含以下信息：

```
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2
chr1	100	.	A	T	60	PASS	.	GT	0/0	0/1
chr1	200	.	G	C	55	PASS	.	GT	0/1	1/1
chr2	300	.	T	A	65	PASS	.	GT	1/1	0/0
```

### 目标ID文件格式

纯文本文件，每行一个变异ID：

```
chr1_100_A_T
chr1_200_G_C
chr2_300_T_A
```

变异ID格式: `{染色体}_{位置}_{参考等位基因}_{替代等位基因}`

## 输出文件格式

### 基因型表 (.gt.txt)

```
CHROM	POS	REF	ALT	sample1	sample2
chr1	100	A	T	AA	AT
chr1	200	G	C	GC	CC
chr2	300	T	A	AA	TT
```

### 序列表 (.seq.txt)

```
CHROM	POS	REF	ALT	sample1	sample2
chr1	100	A	T	A	T
chr1	200	G	C	C	C
chr2	300	T	A	A	T
```

## 性能优化建议

### 内存使用优化

1. **调整批处理大小**:
   - 小文件 (<10MB): `--batch-size 5000`
   - 中等文件 (10-100MB): `--batch-size 10000`
   - 大文件 (>100MB): `--batch-size 50000`

2. **使用压缩输出**:
   ```bash
   python vcf_processor_standalone.py input.vcf targets.txt output --compress
   ```

### 处理速度优化

1. **使用适当的批处理大小**:
   ```bash
   # 对于大文件，增加批处理大小
   python vcf_processor_standalone.py large.vcf targets.txt output --batch-size 50000
   ```

2. **减少日志输出**:
   ```bash
   # 使用静默模式提高处理速度
   python vcf_processor_standalone.py input.vcf targets.txt output --quiet
   ```

## 错误处理

### 常见错误及解决方案

1. **文件不存在错误**
   ```
   错误: VCF文件不存在: input.vcf
   解决: 检查文件路径是否正确
   ```

2. **格式错误**
   ```
   错误: 未找到样本信息
   解决: 检查VCF文件是否包含正确的头部信息
   ```

3. **内存不足**
   ```
   解决: 减少批处理大小 --batch-size 5000
   ```

4. **权限错误**
   ```
   错误: 写入输出文件失败
   解决: 检查输出目录的写入权限
   ```

## 脚本特性

### 支持的VCF格式

- ✅ 标准VCF v4.2格式
- ✅ 压缩VCF文件 (.vcf.gz)
- ✅ 多样本VCF文件
- ✅ 各种基因型格式 (0/0, 0/1, 1/1, ./., 0|1等)

### 基因型转换规则

| 输入基因型 | 参考=A, 替代=T | 输出基因型 | 输出序列 |
|------------|----------------|------------|----------|
| 0/0        | A,A            | AA         | AA       |
| 0/1        | A,T            | AT         | AT       |
| 1/1        | T,T            | TT         | TT       |
| ./.        | 缺失           | NN         | NN       |

### 自定义格式选项

```bash
# 使用VCF标准缺失格式
python vcf_processor_standalone.py input.vcf targets.txt output --miss-fmt "./."

# 使用斜杠分隔基因型
python vcf_processor_standalone.py input.vcf targets.txt output --gt-sep "/"

# 结果: AT 变成 A/T
```

## 与完整版本的对比

| 特性 | 独立脚本版 | 模块化版本 |
|------|------------|------------|
| **命令格式** | `python vcf_processor_standalone.py` | `python -m chip.vcf_processor.cli process` |
| **依赖要求** | 仅Python标准库 | 需要cyvcf2、typer等包 |
| **安装复杂度** | 无需安装 | 需要pip install |
| **处理速度** | 中等 (纯Python) | 更快 (C扩展) |
| **内存效率** | 良好 | 更优 |
| **功能完整性** | 核心功能 | 全部功能 |
| **配置文件支持** | 无 | 支持JSON配置 |
| **多线程支持** | 无 | 支持 |
| **错误处理** | 基础 | 全面 |
| **适用场景** | 快速处理、简单环境 | 生产环境、大规模处理 |

## 选择建议

### 使用独立脚本版本的情况:
- 🎯 **快速原型开发**: 无需复杂环境配置
- 🎯 **教学演示**: 代码简单易懂
- 🎯 **临时处理**: 一次性数据处理任务
- 🎯 **受限环境**: 无法安装额外依赖的环境
- 🎯 **脚本集成**: 容易集成到现有工作流程

### 使用模块化版本的情况:
- 🚀 **生产环境**: 需要高性能和稳定性
- 🚀 **大规模处理**: 处理大型VCF文件
- 🚀 **复杂配置**: 需要配置文件和高级选项
- 🚀 **多线程处理**: 需要并行处理能力
- 🚀 **长期维护**: 项目需要长期维护和扩展

## 总结

`vcf_processor_standalone.py` 是 `python -m chip.vcf_processor.cli process` 的轻量级替代方案，提供相同的核心功能但无需安装依赖。

### 快速选择指南:
- **需要快速开始？** → 使用独立脚本版本
- **需要高性能？** → 使用模块化版本  
- **环境受限？** → 使用独立脚本版本
- **生产环境？** → 使用模块化版本

两种方式的命令参数完全兼容，可以根据需要随时切换！