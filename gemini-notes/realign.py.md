# 代码结构优化-realign

本篇分析 `panel/realign.py` 脚本，该脚本用于探针重比对和注释调整。这是一个典型的"大脚本"，包含 CLI 定义、文件解析、CIGAR 串处理和业务逻辑，非常适合进行模块化重构。

## 1. 核心逻辑抽取：领域模型（Domain Model）
### 📌 问题（Line 31–126）
脚本中包含大量处理 CIGAR 字符串、计算偏移量和提取等位基因的独立函数。这些函数实际上都是围绕"比对记录"这一概念进行的。

```python
def extract_alleles(variant_string): ...
def get_cigar_list(cigar): ...
def get_del_ins(row): ...
def get_pos(row): ...
```

### ✅ 改进：定义 Alignment 类
将这些散落的函数封装到一个 `Alignment` 或 `ProbeMatch` 类中。

```python
@dataclass
class ProbeMatch:
    id: str
    cigar: str
    offset_start: int
    offset_end: int
    strand: str
    # ... 其他字段

    @property
    def cigar_operations(self) -> List[Tuple[int, str]]:
        return re.findall(r"(\d+)([MIDNSHPX=])", self.cigar)

    def calculate_real_position(self, probe_start: int) -> int:
        # 封装 get_pos 的逻辑
        pass

    def get_indel_counts(self) -> InsDelCount:
        # 封装 get_del_ins 的逻辑
        pass
```

## 2. 外部工具封装：适配器模式（Adapter Pattern）
### 📌 问题（Line 206–231）
脚本直接使用 `delegator.run` 调用 `bedtools` 和 `minimap2`，命令字符串拼接散落在各个函数中。

```python
cmd_line = f"bedtools slop -b {flank_size} -i {target_bed} -g {genome_fai} > {flank_bed}"
delegator.run(cmd_line)
```

### ✅ 改进：工具包装器
创建专门的类来封装外部工具调用，提供类型安全的方法接口。

```python
class BedTools:
    def slop(self, input_bed: Path, genome_file: Path, flank: int, output: Path):
        cmd = ["bedtools", "slop", "-b", str(flank), "-i", str(input_bed), "-g", str(genome_file)]
        # 执行命令并处理错误
        
class Minimap2:
    def align(self, target: Path, query: Path, output: Path, threads: int = 4):
        # ...
```

## 3. CLI 组织：命令组与子命令
### 📌 问题（Line 300–605）
`realign`, `realign2`, `realign3` 等多个命令逻辑高度相似，存在大量重复代码（如生成 flank bed, fasta, paf 的流程）。

### ✅ 改进：Pipeline 模式
将通用的处理流程抽象为 Pipeline，不同的命令只是 Pipeline 的不同配置或入口。

```python
class RealignPipeline:
    def __init__(self, config: RealignConfig):
        self.config = config
    
    def run(self):
        self._prepare_bed()
        self._generate_fasta()
        self._run_minimap()
        self._process_paf()

@app.command()
def vcf_mode(...):
    config = RealignConfig(input_type='vcf', ...)
    RealignPipeline(config).run()

@app.command()
def bed_mode(...):
    config = RealignConfig(input_type='bed', ...)
    RealignPipeline(config).run()
```

## 🎯 改进总结
### 提升点
*   **内聚性**：将数据和操作数据的逻辑封装在一起（`ProbeMatch`），代码更易理解。
*   **可测试性**：外部工具封装后，可以轻松 mock `BedTools` 类进行单元测试，而不需要实际安装工具。
*   **复用性**：Pipeline 模式消除了 `realign` 和 `realign2` 之间的重复代码。

### 适用场景
*   复杂的生物信息学流程脚本。
*   需要频繁调用外部命令行工具的项目。
