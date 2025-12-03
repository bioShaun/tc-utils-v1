# Group Map 逻辑修复说明

## 问题描述

**发现问题：** Group map文件的格式理解有误，导致映射逻辑错误。

**Group Map文件格式：**
- 输入格式：`group_id\tgene_id(ref_gene_id)`
- 即：第二列是 `ref_gene_id`，不是 `qry_gene_id`

## 修复内容

### 1. extract_gene_id_map_from_tmap.py

**修复前（第121行）：**
```python
for qry_gene, ref_gene, _, _ in mappings:
    if qry_gene in gene_to_group:  # ❌ 错误：应该用ref_gene查找
        group_id = gene_to_group[qry_gene]
        f.write(f"{group_id}\t{qry_gene}\n")
```

**修复后（第122行）：**
```python
for qry_gene, ref_gene, _, _ in mappings:
    # group_map是ref_gene_id到group_id的映射，所以用ref_gene查找
    if ref_gene in gene_to_group:  # ✅ 正确
        group_id = gene_to_group[ref_gene]
        f.write(f"{group_id}\t{qry_gene}\n")
```

### 2. filter_gene_map_by_chromosome.py

**修复前（第149-221行）：**
- 使用 `valid_qry_genes` 集合来过滤group_map
- 错误地认为group_map第二列是qry_gene
- 没有建立ref_gene到qry_gene的映射关系

**修复后（第149-221行）：**
```python
def filter_group_map(group_map_file, gene_map_file, qry_to_chrom, ref_to_chrom, genome_mapping, gene_to_group):
    """
    注意：group_map文件格式是 group_id -> ref_gene_id
    """
    # 首先建立ref_gene到qry_gene的映射，并过滤染色体不匹配的记录
    ref_to_qry_valid = {}

    with open(gene_map_file, "r") as f:
        header = f.readline().strip()

        for line in f:
            # ... 染色体检查逻辑 ...
            # 通过所有检查，建立ref_gene到qry_gene的映射
            ref_to_qry_valid[ref_gene] = qry_gene

    # 过滤group map
    with open(group_map_file, "r") as f:
        header = f.readline().strip()

        for line in f:
            fields = line.strip().split("\t")

            group_id = fields[0]
            ref_gene = fields[1]  # group_map中的gene_id实际上是ref_gene_id

            # 检查这个ref_gene是否对应有效的qry_gene
            if ref_gene in ref_to_qry_valid:
                qry_gene = ref_to_qry_valid[ref_gene]
                filtered.append((group_id, qry_gene))
```

### 3. 测试数据修复

**所有测试文件中的GROUP_MAP_CONTENT都已修正：**

```python
# 修复前（错误）
GROUP_MAP_CONTENT = """group_id\tgene_id
1\tTraesCS1A01G000100  # 错误：应该是ref_gene
2\tTraesCS1A01G000100LC  # 错误：应该是ref_gene
"""

# 修复后（正确）
GROUP_MAP_CONTENT = """group_id\tgene_id
1\tTraesCS1A01G000100  # 正确：这是ref_gene
2\tTraesCS1A01G000200  # 正确：这是ref_gene
"""
```

## 逻辑说明

### 正确的映射流程

1. **tmap文件**：
   ```
   ref_gene_id         qry_gene_id
   TraesCS1A01G000100  TraesCS1A01G000100
   TraesCS1A01G000100  TraesCS1A01G000100LC
   ```

2. **group_map文件**：
   ```
   group_id  gene_id(ref_gene_id)
   1         TraesCS1A01G000100
   ```

3. **输出文件**：
   ```
   group_id  qry_gene_id
   1         TraesCS1A01G000100
   1         TraesCS1A01G000100LC
   ```

### 为什么这样设计？

- group_map文件描述的是 ref_gene 到 group 的关系
- 一个 ref_gene 可能对应多个 qry_gene（不同转录本、异构体等）
- 输出时需要输出每个 qry_gene 对应的 group_id

## 测试验证

✅ **所有测试通过：**
- 原生Python测试：7/7 通过
- pytest测试：7/7 通过

### 测试覆盖率

1. **extract_gene_id_map_from_tmap.py**：
   - ✅ 不使用group map
   - ✅ 使用group map（修正后）
   - ✅ 优先级处理

2. **filter_gene_map_by_chromosome.py**：
   - ✅ 不使用group map的过滤
   - ✅ 使用group map的过滤（修正后）
   - ✅ 染色体不匹配的情况
   - ✅ 基因在GFF中不存在的情况

## 影响范围

- ✅ 不影响现有功能（向后兼容）
- ✅ 只修正了group_map的查找逻辑
- ✅ 测试数据已同步修正
- ✅ 文档已更新

## 文件修改列表

1. `extract_gene_id_map_from_tmap.py` - 修正group_map查找逻辑
2. `filter_gene_map_by_chromosome.py` - 修正filter_group_map函数
3. `test_extract_gene_id_map.py` - 修正测试数据
4. `test_extract_gene_id_map_pytest.py` - 修正测试数据
5. `test_filter_gene_map_by_chromosome.py` - 修正测试数据
6. `test_filter_gene_map_by_chromosome_pytest.py` - 修正测试数据

## 修复时间

2025年12月3日

## 总结

此次修复澄清了group_map文件的格式，并修正了两个脚本中的查找逻辑。修复后的代码能够正确处理group_map文件，实现从ref_gene到group_id再到qry_gene的正确映射关系。
