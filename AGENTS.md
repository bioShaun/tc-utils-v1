# tc-pytools

生物信息学数据处理（GTF/BED/VCF 文件操作）和系统自动化的 Python CLI 工具集合。

## 项目结构
- `gtf/` - 基因组文件操作工具
  - `split_bed.py` - 旧版脚本（仅供参考，**禁止**作为模板）
  - `split_large_genome_bed_v2.py` - **标准模板**（新脚本必须参考此文件！）
- `tests/` - Pytest 测试套件
- `pyproject.toml` - 依赖管理

## 可用技能

需要时使用 `skill` 工具加载详细说明：

| 技能 | 触发场景 | 用法 |
|------|----------|------|
| `python-review` | 代码审查、安全分析、质量评估 | `skill(name="python-review")` |
| `python-refactor` | 重构遗留代码、应用整洁代码原则 | `skill(name="python-refactor")` |
| `python-modern-cli` | 从零创建新的 CLI 工具 | `skill(name="python-modern-cli")` |

## 外部文件加载

**重要**：处理此项目时，使用延迟加载参考文件：
- 重构或创建 CLI 工具时，读取 `@gtf/split_large_genome_bed_v2.py` 作为参考实现
- **不要**预先加载所有文件——仅在特定任务需要时加载

## 代码规范

### 必需技术栈
- **CLI 框架：** `typer` 配合 `Annotated` 类型提示
- **日志：** `loguru`（**禁止**使用 `print` 或标准 `logging`）
- **终端输出：** `rich.console.Console` 用于用户交互信息
- **数据处理：** `polars`（首选）或 `pandas`（仅限遗留/小数据）
- **路径：** `pathlib.Path`（**禁止**使用 `os.path`）
- **生物信息：** `cyvcf2`（VCF）、`pyfaidx`（FASTA）、`pysam`（BAM/SAM）

### 文档规范

**所有脚本必须包含：**

1. **中文模块文档字符串**：脚本开头说明功能、使用示例和输出格式
2. **中文函数文档字符串**：所有公开函数使用 Google/NumPy 风格
3. **中文注释**：关键逻辑使用中文注释说明

**模块文档字符串模板：**
```python
"""
脚本功能简述。

详细说明脚本的用途和特点。

使用示例:
    # 基本用法
    python script.py input.vcf output.csv

    # 带可选参数
    python script.py input.vcf output.csv -c config.txt -v

输出格式:
    说明输出文件的格式和各列含义。
"""
```

**函数文档字符串模板：**
```python
def process_data(input_path: Path, threshold: float = 0.5) -> pl.DataFrame:
    """
    处理输入数据并返回结果。

    Args:
        input_path: 输入文件路径
        threshold: 过滤阈值，默认 0.5

    Returns:
        处理后的 DataFrame

    Raises:
        typer.Exit: 输入文件不存在时退出
    """
```

### 生物信息学库

使用现代、活跃维护的包处理生物信息学数据：

| 数据类型 | 首选包 | 使用场景 | 备选方案 |
|----------|--------|----------|----------|
| **VCF 文件** | `cyvcf2` | 所有 VCF 解析/读取任务 | `pysam.VariantFile`（cyvcf2 缺少特定功能时） |
| **FASTA 文件** | `pyfaidx` | 大基因组、随机访问 | `Biopython SeqIO`（小文件、格式转换、序列操作） |
| **BAM/SAM 文件** | `pysam` | 所有比对文件操作 | 无 |
| **BED/GTF 文件** | `polars` / `pandas` | 表格型基因组坐标 | 自定义解析器（尽量避免） |

**版本要求：** 见 `pyproject.toml` 中的最低支持版本。

**VCF 处理：**
```python
# 首选：cyvcf2 高性能读取
from cyvcf2 import VCF

for variant in VCF('variants.vcf.gz'):
    chrom, pos = variant.CHROM, variant.POS
    ref, alt = variant.REF, variant.ALT[0]
    
    # 高效访问基因型
    gts = variant.gt_types  # 0=HOM_REF, 1=HET, 2=HOM_ALT, 3=UNKNOWN
    
# 备选：pysam.VariantFile（cyvcf2 功能不足时）
from pysam import VariantFile

with VariantFile('output.vcf', 'w', header=custom_header) as vcf_out:
    vcf_out.write(record)
```

**FASTA 处理：**
```python
# 首选：pyfaidx 处理大基因组（内存高效、索引访问）
from pyfaidx import Fasta

genome = Fasta('genome.fa')
sequence = genome['chr1'][1000:2000]  # O(1) 随机访问
reverse_comp = genome['chr1'][1000:2000].reverse.complement

# 可接受：Biopython 特定场景
from Bio import SeqIO
from Bio.Seq import Seq

# ✓ 适用：小文件、格式转换
records = list(SeqIO.parse('small.fasta', 'fasta'))
SeqIO.convert('input.fasta', 'fasta', 'output.gb', 'genbank')

# ✓ 适用：序列操作
protein = Seq('ATGGCCATTGTAATG').translate()

# ✗ 避免：大基因组解析（改用 pyfaidx）
genome = {rec.id: rec.seq for rec in SeqIO.parse('genome.fa', 'fasta')}  # 内存效率低！
```

**BAM/SAM 处理：**
```python
import pysam

# 区域提取（内存高效）
with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for read in bam.fetch('chr1', 1000, 2000):
        if read.mapping_quality >= 30:
            print(f"{read.query_name}: {read.reference_start}")
            
# 索引操作
pysam.index('input.bam')
```

**性能指南：**
- **cyvcf2** 处理大 VCF 文件比 PyVCF 快 5-10 倍
- **pyfaidx** 无论基因组大小都占用极少内存（通过 FAIDX 索引）
- **pysam** 高效处理 BAM/CRAM，内置解压缩

**Biopython 适用场景：**
- 文件大小 < 100MB
- 格式转换（FASTA ↔ GenBank 等）
- 序列操作（翻译、反向互补、motif 查找）
- 系统发育和比对任务（Phylo、AlignIO 模块）

### 类型提示（必需）
```python
# 使用 Python 3.10+ 语法
def process(items: list[str], config: dict[str, Any] | None = None) -> Path:
    ...
```

### 文档字符串（必需）
所有公开函数使用 Google 或 NumPy 风格：
```python
def validate_file(path: Path) -> bool:
    """
    验证文件是否存在且可读。
    
    Args:
        path: 要验证的文件路径
        
    Returns:
        文件有效返回 True
        
    Raises:
        typer.Exit: 验证失败时退出
    """
```

## 错误处理模式

所有 CLI 工具遵循此标准模式：

```python
from loguru import logger
from rich.console import Console

console = Console()

def validate_input(path: Path) -> None:
    """验证输入，正确处理错误。"""
    if not path.exists():
        logger.error(f"文件不存在: {path}")
        console.print(f"[red]错误:[/red] 文件不存在: {path}")
        raise typer.Exit(code=1)
```

**规则：**
- 使用 `logger.error()` 记录错误用于调试
- 使用 `console.print()` 配合 Rich 标记显示用户友好信息
- 错误使用 `[red]错误:[/red]`，警告使用 `[yellow]警告:[/yellow]`，成功使用 `[green]✓[/green]`
- 失败时始终 `raise typer.Exit(code=1)`，**禁止** `sys.exit()`

## 测试规范

### 框架和命令
```bash
pytest tests/                    # 运行所有测试
pytest tests/ -v                 # 详细输出
pytest tests/ --cov=gtf          # 带覆盖率
```

### 命名约定
- 测试文件：`test_<模块名>.py`
- 测试函数：`test_<函数名>_<场景>()`
- 夹具：使用 `conftest.py` 存放共享夹具

### 结构
```python
def test_validate_bed_file_returns_error_for_missing_file(tmp_path: Path) -> None:
    """测试对不存在文件的验证失败。"""
    # 准备
    fake_path = tmp_path / "nonexistent.bed"
    
    # 执行和断言
    with pytest.raises(SystemExit):
        validate_bed_file(fake_path)
```

## Git 约定

### 提交信息格式
```
<类型>: <简短描述>

<可选正文>
```

**类型：**
- `feat`: 新功能
- `fix`: 修复 bug
- `refactor`: 代码重构（无功能变化）
- `docs`: 仅文档
- `test`: 添加或更新测试
- `chore`: 维护任务

**示例：**
```
feat: 添加分割 BED 文件的坐标转换功能
refactor: 使用 typer 和 loguru 现代化 split_bed.py
fix: 优雅处理空 BED 文件
```

## 开发工作流

### 重构遗留脚本
遇到使用 `argparse`、`os.path` 或 `print` 的脚本时：
1. 加载 `python-refactor` 技能
2. 读取 `gtf/split_large_genome_bed_v2.py` 作为参考
3. **分析**现有逻辑，识别反模式
4. 使用 `typer`、`loguru`、`rich` 和 `pathlib` **现代化**
5. 显式**验证**所有输入，提供正确的错误信息
6. **测试**重构后的脚本

### 创建新工具
1. 加载 `python-modern-cli` 技能
2. 读取 `gtf/split_large_genome_bed_v2.py` 作为参考
3. 使用 `Annotated` 参数创建 `typer` 骨架
4. 尽早添加验证函数
5. 在 `main()` 中根据 `--verbose` 标志设置日志
6. 同步编写测试

## 常见问题

| 问题 | 解决方案 |
|------|----------|
| `typer`/`pandas` 导入错误 | 检查 `pyproject.toml` 依赖 |
| 处理 >1GB 文件很慢 | 从 `pandas` 切换到 `polars` |
| 类型提示不工作 | 确保 Python 3.10+ 和正确导入 |
| Rich 颜色不显示 | 检查终端是否支持 ANSI 颜色 |
| Biopython SeqIO 处理大 FASTA 很慢 | 使用 `pyfaidx` 进行索引随机访问 |
| VCF 解析内存错误 | 切换到 `cyvcf2` 配合流式迭代 |
| BAM 文件崩溃或卡住 | 使用 `pysam` 的区域 `fetch()` 而非迭代所有读段 |
