# compareGT_v3.py 优化总结

## 概述

将 VCF 基因型比较工具从 pandas 字符串处理升级为 numpy 向量化矩阵运算，实现 20+ 倍性能提升。

## 性能提升

| 指标 | V2 (pandas) | V3 (numpy) | 提升 |
|------|-------------|------------|------|
| 300 对比较 | 4分9秒 | **11秒** | **~22x** |
| 62481 对比较 | >10分钟 | **1分36秒** | **>6x** |
| 内存占用 | ~2-4 GB | ~3.5 MB/块 | **~100x** |

## 关键优化策略

### 1. 数据结构优化

```python
# 之前：字符串处理（慢，内存大）
gt = variant.gt_bases[i]  # "A/G", "./."
gt = gt.replace("|", "/")

# 之后：数值编码（快，内存小）
gt = variant.gt_types  # numpy int8 数组
# 0=HOM_REF, 1=HET, 2=UNKNOWN, 3=HOM_ALT
```

- `gt_bases` 字符串 → `gt_types` int8 数值编码
- polars DataFrame → numpy 矩阵
- 内存占用从 GB 级降到 MB 级

### 2. 算法优化

```python
# 之前：Python 循环逐对计算
for pair in pairs:  # 62481 次循环/块
    stats = compute_chunk_stats(chunk_df, pair)

# 之后：numpy 广播一次性计算所有对
gt_a = gt_matrix[:, idx_a].T  # (n_pairs, n_variants)
gt_b = gt_matrix[:, idx_b].T
gt_equal = (gt_a == gt_b) & valid  # 向量化比较
```

- Python 循环 → numpy 广播向量化
- 逐对计算 → 所有样本对同时计算
- 重复分类 → 一次性预计算

### 3. 分块处理

```python
def load_vcf_as_matrix(vcf_path, chunk_size=10000):
    """分块读取，内存可控。"""
    for variant in vcf:
        chunk.append(variant.gt_types.copy())
        if len(chunk) >= chunk_size:
            yield np.array(chunk, dtype=np.int8)
            chunk = []
```

- 支持任意大小 VCF（内存可控）
- 累加器模式聚合统计量
- `--chunk-size` 参数可调

## 踩坑经验

### 1. gt_types 编码

cyvcf2 的 `gt_types` 编码：

| 值 | 含义 | 基因型示例 |
|---|------|-----------|
| 0 | HOM_REF | 0/0 |
| 1 | HET | 0/1, 1/0 |
| 2 | UNKNOWN | ./. |
| 3 | HOM_ALT | 1/1 |

**注意**：不是直觉的 `0,1,2,3=HOM_REF,HET,HOM_ALT,UNKNOWN`

### 2. Multi-allelic 位点处理

**问题**：`gt_types` 只区分 HOM/HET，不区分不同杂合等位基因

```
# Multi-allelic 位点 (ALT=G,C)
Sample1: A/G (0/1) -> gt_types=1 (HET)
Sample2: A/C (0/2) -> gt_types=1 (HET)
# gt_types 相同，但实际基因型不同！
```

**解决方案**：将 multi-allelic 位点标记为 UNKNOWN，排除出有效位点计算

```python
if len(variant.ALT) > 1:
    gt[:] = GT_UNKNOWN  # 标记为缺失
```

- 仅占总位点的 ~0.36%，影响可忽略
- 总位点数保持不变，有效位点排除 multi-allelic

### 3. 结果一致性验证

排除 multi-allelic 后，V2 和 V3 的关键统计量完全一致：

```
V2 有效位点=115855, 相似度%=86.954
V3 有效位点=115855, 相似度%=86.954
```

## 代码规范

- 中文注释和文档字符串
- `Annotated` 类型提示
- Rich 进度条 + loguru 日志
- 完整的单元测试覆盖（22 个测试）

## 使用示例

```bash
# 比较所有样本对
python compareGT_v3.py input.vcf output.csv

# 指定比较列表
python compareGT_v3.py input.vcf output.csv -c compare_pairs.txt

# 调整分块大小（内存不足时减小）
python compareGT_v3.py input.vcf output.csv --chunk-size 5000

# 详细日志
python compareGT_v3.py input.vcf output.csv -v
```

## 总结

| 优化点 | 效果 |
|--------|------|
| gt_types 数值编码 | 避免字符串解析 |
| numpy 向量化 | 消除 Python 循环 |
| 分块处理 | 内存可控 |
| 排除 multi-allelic | 保证结果正确性 |

**核心经验**：对于大规模数据处理，优先使用数值编码 + numpy 向量化，避免 Python 层面的循环和字符串操作。
