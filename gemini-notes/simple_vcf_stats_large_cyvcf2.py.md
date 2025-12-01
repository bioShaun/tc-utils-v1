# 代码结构优化-simple_vcf_stats_large_cyvcf2

本篇分析 `panel/simple_vcf_stats_large_cyvcf2.py` 脚本，该脚本用于并行处理 VCF 文件。脚本手动管理了进程池、批处理和 Future 对象，导致并发控制逻辑非常复杂且难以维护。

## 1. 并发控制：简化并行模式
### 📌 问题（Line 253–321）
`process_variants_in_parallel` 函数手动实现了复杂的生产者-消费者模式：手动分批、提交任务、管理 `active_futures`、轮询完成状态。这不仅代码量大（近 70 行），而且容易出现死锁或资源泄露。

```python
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        # ... 初始化变量
        for variant in vcf:
            # ... 手动分批
            if len(batch) >= batch_size:
                # ... 手动流控
                while active_futures >= max_active_futures:
                    # ... 轮询检查
```

### ✅ 改进：使用高级并发工具
利用 `concurrent.futures.as_completed` 或 `tqdm` 的并发包装器来简化逻辑。Python 的迭代器机制可以自动处理流式数据。

```python
from itertools import islice

def chunked_iterable(iterable, size):
    it = iter(iterable)
    while True:
        chunk = list(islice(it, size))
        if not chunk:
            break
        yield chunk

def process_variants_in_parallel_simplified(vcf, writer, num_processes, batch_size):
    # 将 VCF 迭代器转换为数据块生成器
    # 注意：这里需要处理 cyvcf2 对象的序列化问题，通常在主进程转为 dict/dataclass
    data_stream = (variant_to_data(v) for v in vcf)
    batches = chunked_iterable(data_stream, batch_size)

    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        # 提交所有任务（注意：对于无限流，需控制提交速度，这里假设文件有限或使用 map）
        # 使用 map 可以保持顺序，且代码最简洁
        results_iterator = executor.map(process_variant_batch, batches)
        
        for batch_results in results_iterator:
            write_results(writer, batch_results)
```

## 2. 数据传输：减少序列化开销
### 📌 问题（Line 58–76）
为了在进程间传递数据，脚本将 `cyvcf2.Variant`（C 扩展对象）转换为 Python `VariantData` 对象。这涉及大量属性访问和列表创建，对于数百万行变异位点，序列化/反序列化（Pickle）开销巨大。

```python
def variant_to_data(variant):
    return VariantData(
        chrom=variant.CHROM,
        gt_types=variant.gt_types.tolist(), # numpy 转 list
        # ...
    )
```

### ✅ 改进：使用共享内存或轻量级结构
对于数值型数据（如基因型），尽量保持 Numpy 数组形态，或者仅传递必要字段的元组（Tuple），比 Dataclass 更快。或者考虑使用 `joblib` 的内存映射功能。

```python
# 仅提取计算所需的最小字段
def variant_to_tuple(variant):
    return (variant.CHROM, variant.POS, variant.gt_types, variant.aaf)

# 在子进程中重建上下文或直接计算
def process_batch_tuples(batch_tuples):
    results = []
    for chrom, pos, gt_types, aaf in batch_tuples:
        # 直接计算，无需重建完整对象
        pass
    return results
```

## 3. 进度监控：解耦 UI 逻辑
### 📌 问题（Line 189–200）
进度打印逻辑硬编码在处理循环中，且直接操作 `sys.stdout`。

### ✅ 改进：使用 tqdm
`tqdm` 是 Python 标准的进度条库，能自动处理平滑更新、预计剩余时间等。

```python
from tqdm import tqdm

# 配合 executor.map 使用
results = list(tqdm(executor.map(func, data), total=total_count))
```

## 🎯 改进总结
### 提升点
*   **简洁性**：并发核心逻辑从 70 行减少到 10 行左右。
*   **性能**：减少不必要的对象转换和序列化开销。
*   **可读性**：使用标准库惯用法（`executor.map`），降低认知负担。

### 适用场景
*   CPU 密集型的大文件处理任务。
*   需要高性能并行计算的数据清洗脚本。
