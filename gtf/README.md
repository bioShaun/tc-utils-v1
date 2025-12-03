# 基因映射工具

本目录包含两个用于处理基因映射的Python脚本及相关测试。

## 文件结构

```
gtf/
├── extract_gene_id_map_from_tmap.py      # 基因映射提取脚本
├── filter_gene_map_by_chromosome.py      # 染色体过滤脚本
├── README.md                             # 文档说明
└── tests/                                # 测试目录
    ├── __init__.py
    ├── test_extract_gene_id_map.py       # 测试extract_gene_id_map_from_tmap.py
    ├── test_filter_gene_map_by_chromosome.py  # 测试filter_gene_map_by_chromosome.py
    └── run_all_tests.py                   # 运行所有测试
```

## 脚本详情

### 1. extract_gene_id_map_from_tmap.py

从gffcompare的.tmap文件中提取基因ID映射。

**功能:**
- 解析tmap文件，提取qry_gene_id与ref_gene_id的映射关系
- 根据class code优先级过滤（= > c ≈ k > m > j）
- 可选：根据group map文件生成group_id与qry_gene_id的映射

**使用示例:**
```bash
# 基本使用
python extract_gene_id_map_from_tmap.py input.tmap -o output.txt

# 使用group map
python extract_gene_id_map_from_tmap.py input.tmap -g group_map.txt -o output.txt
```

**输出文件:**
- `gene_id_map.txt`: qry_gene_id\tref_gene_id\tclass_code
- `gene_id_map_group_map.txt`: group_id\tgene_id (仅当使用-g参数时)

---

### 2. filter_gene_map_by_chromosome.py

根据染色体信息过滤基因映射。

**功能:**
- 从GFF文件中提取基因ID到染色体ID的映射
- 读取genome映射文件（qry_genome,ref_genome）
- 过滤基因映射，只保留染色体匹配的记录
- 可选：同时过滤group map文件

**使用示例:**
```bash
# 过滤基因映射
python filter_gene_map_by_chromosome.py \
    --gene-map gene_id_map.txt \
    --qry-gff query_genes.gff \
    --ref-gff ref_genes.gff \
    --genome-map genome_map.txt \
    --output filtered.txt

# 同时过滤group map
python filter_gene_map_by_chromosome.py \
    --gene-map gene_id_map.txt \
    --group-map gene_id_map_group_map.txt \
    --qry-gff query_genes.gff \
    --ref-gff ref_genes.gff \
    --genome-map genome_map.txt \
    --output filtered.txt
```

**输出文件:**
- `filtered_gene_map.txt`: 过滤后的基因映射
- `filtered_gene_map_filtered_group_map.txt`: 过滤后的group映射（仅当使用--group-map时）

---

## 完整工作流程

```bash
# 步骤1：从tmap提取基因映射
python extract_gene_id_map_from_tmap.py input.tmap -g group_map.txt -o gene_id_map.txt

# 步骤2：根据染色体过滤
python filter_gene_map_by_chromosome.py \
    --gene-map gene_id_map.txt \
    --group-map gene_id_map_group_map.txt \
    --qry-gff query_genes.gff \
    --ref-gff ref_genes.gff \
    --genome-map genome_map.txt \
    --output final_filtered.txt

# 最终输出：
# - final_filtered.txt: 过滤后的基因映射
# - final_filtered_filtered_group_map.txt: 过滤后的group映射
```

---

## 文件格式说明

### tmap文件
gffcompare生成的.tmap文件，包含基因对应关系和class code。

### group_map.txt
```
group_id\tgene_id
1\tTraesCS1A01G000100
2\tTraesCS1A01G000100LC
```

### genome_map.txt
```
qry_genome,ref_genome
chr1,RefChr1
chr2,RefChr2
```

### GFF/GTF文件
标准GFF/GTF格式，脚本支持多种基因ID提取格式：
- `gene_id="xxxxx"`
- `gene_id = "xxxxx"`
- `gene_id "xxxxx"`
- `ID=xxxxx`
- `Name=xxxxx`

---

## 测试

测试文件位于`tests/`目录下：

```bash
cd tests

# 运行单个测试
python test_extract_gene_id_map.py
python test_filter_gene_map_by_chromosome.py

# 或运行所有测试
python run_all_tests.py
```

测试覆盖：
- 基本功能
- 边界情况
- 错误处理
- 文件格式兼容性

测试结果：7个测试用例全部通过 ✓

---

## 版本历史

- v1.0: 初始版本，支持tmap文件解析和基因映射提取
- v1.1: 增加group map支持和染色体过滤功能
- v1.2: 增加完整测试套件
- v1.3: 重构测试结构，将测试文件移至tests/目录
