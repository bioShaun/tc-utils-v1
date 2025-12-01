# tsv2excel.py - 通用文件转换工具分析

## 代码概览

**文件路径**: utils/tsv2excel.py

**功能描述**: 将CSV或TSV文件转换为Excel格式，支持自动检测分隔符。

**核心价值**: 简单的工具类函数，展示了**最小可行产品(MVP)**的设计思路。

## 代码结构分析

### 设计特点

**极简主义设计**:
- 单一函数完成核心功能
- 最少的依赖项
- 直接的逻辑流程

### 流程图

```
输入文件验证 → 格式检测 → 读取数据 → 保存Excel → 返回结果
```

### 优点

1. **清晰的函数签名**
   ```python
   def file2excel(input_path: Path, excel_path: Path) -> Tuple[int, int]:
   ```
   - 明确的输入输出
   - 完整的类型注解
   - 语义化的参数名

2. **健壮的格式检测** (第29-50行)
   ```python
   # 方法1: 基于扩展名
   if input_path.suffix.lower() == '.tsv':
       df = pd.read_table(input_path)
   elif input_path.suffix.lower() == '.csv':
       df = pd.read_csv(input_path)

   # 方法2: 基于内容(扩展名无法识别时)
   else:
       with open(input_path, 'r', encoding='utf-8') as f:
           first_line = f.readline()

       if '\t' in first_line:
           df = pd.read_table(input_path)
       elif ',' in first_line:
           df = pd.read_csv(input_path)

   # 方法3: Python引擎自动检测
   else:
       df = pd.read_csv(input_path, sep=None, engine='python')
   ```
   - 三层检测策略
   - 渐进式回退
   - 最坏情况下使用自动检测

3. **完整的错误处理**
   - 文件存在性检查
   - 详细的异常消息
   - 异常链传递(`raise ... from e`)

4. **用户友好的日志**
   ```python
   logger.info(f"Reading file: {input_path}")
   logger.info(f"Data loaded. Shape: {df.shape}")
   logger.info(f"Successfully saved to {excel_path}")
   ```
   - 关键步骤的日志记录
   - 包含上下文信息

5. **返回值设计**
   ```python
   return df.shape
   ```
   - 返回数据维度，便于调用方验证
   - 元组形式直观易用

### 技术细节

1. **文件编码处理**
   ```python
   with open(input_path, 'r', encoding='utf-8') as f:
       first_line = f.readline()
   ```
   - 明确指定UTF-8编码
   - 读取首行进行格式检测

2. **异常链传递**
   ```python
   raise ValueError(f"Could not determine file format for {input_path}") from e
   ```
   - 使用`from e`保留原始异常信息
   - 便于调试和问题追踪

3. **pandas API使用**
   - `read_table()`: TSV文件
   - `read_csv()`: CSV文件
   - `to_excel()`: 输出Excel文件
   - `index=False`: 不保存索引列

## 可维护性分析

### 优秀实践

1. **函数单一职责**
   - 只做一件事：文件格式转换
   - 易于测试和复用

2. **清晰的文档字符串**
   ```python
   """
   Convert a CSV or TSV file to an Excel file.

   Args:
       input_path (Path): Path to the input CSV or TSV file.
       excel_path (Path): Path to the output Excel file.

   Returns:
       Tuple[int, int]: The shape of the dataframe (rows, columns).

   Raises:
       FileNotFoundError: If the input file does not exist.
       ValueError: If the file format cannot be determined or parsed.
   """
   ```
   - 完整的参数说明
   - 返回值说明
   - 异常说明
   - 符合Google风格

3. **类型注解完整**
   - 所有参数和返回值都有类型提示
   - 使用`Path`和`Tuple`等语义化类型

### 潜在问题

1. **硬编码编码格式**
   ```python
   encoding='utf-8'
   ```
   - **影响**: 无法处理非UTF-8编码文件
   - **建议**: 支持编码参数或自动检测

2. **pandas依赖**
   ```python
   import pandas as pd
   ```
   - **影响**: 对于简单转换可能过于重量级
   - **建议**: 考虑使用`csv`模块作为轻量级替代

3. **Excel格式固定**
   ```python
   df.to_excel(excel_path, index=False)
   ```
   - **影响**: 无法自定义工作表名、格式等
   - **建议**: 支持额外参数配置

4. **内存使用**
   - 将整个文件加载到内存
   - **影响**: 大文件处理可能内存不足
   - **建议**: 支持分块读取

### 扩展性限制

1. **输出格式单一**
   - 仅支持Excel(.xlsx)
   - 无法输出其他格式(CSV、TSV、JSON等)

2. **转换规则固定**
   - 无法自定义转换逻辑
   - 无法添加计算列

3. **批量转换不支持**
   - 只能处理单个文件
   - 无法批量转换多个文件

## 改进方案

### 1. 配置增强

```python
@dataclass
class ConvertConfig:
    """转换配置"""
    input_encoding: str = "utf-8"
    output_format: str = "xlsx"  # xlsx, xls, csv, tsv
    include_index: bool = False
    sheet_name: str = "Sheet1"
    float_precision: int = 3

    @staticmethod
    def auto_detect_encoding(file_path: Path) -> str:
        """自动检测文件编码"""
        # 使用chardet库检测编码
        import chardet
        with open(file_path, 'rb') as f:
            raw_data = f.read(10000)  # 读取前10KB
            result = chardet.detect(raw_data)
            return result['encoding'] or 'utf-8'
```

### 2. 批量转换

```python
def batch_convert(
    input_dir: Path,
    output_dir: Path,
    pattern: str = "*.{csv,tsv}",
    recursive: bool = False
) -> List[Tuple[Path, Path, Tuple[int, int]]]:
    """批量转换目录中的所有文件"""
    results = []
    output_dir.mkdir(parents=True, exist_ok=True)

    for input_file in input_dir.glob(pattern):
        output_file = output_dir / f"{input_file.stem}.xlsx"
        try:
            shape = file2excel(input_file, output_file)
            results.append((input_file, output_file, shape))
            logger.success(f"Converted: {input_file} → {output_file}")
        except Exception as e:
            logger.error(f"Failed to convert {input_file}: {e}")
            results.append((input_file, output_file, (0, 0)))

    return results
```

### 3. 管道模式

```python
from typing import Callable

class DataProcessor:
    """数据处理器，支持管道式转换"""

    def __init__(self, df: pl.DataFrame):
        self.df = df
        self.transformers: List[Callable] = []

    def add_column(self, name: str, expression: str):
        """添加计算列"""
        def transform(df):
            return df.with_columns(pl.col(expression).alias(name))
        self.transformers.append(transform)
        return self

    def filter_rows(self, condition: str):
        """添加行过滤条件"""
        def transform(df):
            return df.filter(pl.col(condition))
        self.transformers.append(transform)
        return self

    def execute(self) -> pl.DataFrame:
        """执行所有转换"""
        for transformer in self.transformers:
            self.df = transformer(self.df)
        return self.df

# 使用示例
processor = DataProcessor(df)
result = (processor
    .add_column("total", "col1 + col2")
    .filter_rows("total > 100")
    .execute())
```

### 4. 流式处理

```python
def convert_large_file(
    input_path: Path,
    output_path: Path,
    chunk_size: int = 10000
) -> Tuple[int, int]:
    """处理大文件的流式转换"""
    total_rows = 0
    total_cols = None

    with pd.ExcelWriter(output_path) as writer:
        for chunk in pd.read_csv(input_path, chunksize=chunk_size):
            if total_cols is None:
                total_cols = chunk.shape[1]

            # 写入Excel，支持多工作表
            sheet_name = f"Chunk_{total_rows // chunk_size + 1}"
            chunk.to_excel(writer, sheet_name=sheet_name, index=False)

            total_rows += len(chunk)
            logger.info(f"Processed {total_rows} rows...")

    return total_rows, total_cols or 0
```

## 学习要点

### Python技巧

1. **异常链**
   ```python
   raise ValueError(...) from e
   ```
   - 保留原始异常信息
   - 便于调试

2. **类型注解**
   ```python
   def func() -> Tuple[int, int]:
   ```
   - 提升代码可读性
   - IDE自动补全

3. **文档字符串**
   - 描述函数功能
   - 说明参数和返回值
   - 列出可能抛出的异常

### 最佳实践

1. **防御性编程**
   - 检查文件存在性
   - 验证输入参数
   - 处理边界情况

2. **用户友好**
   - 清晰的错误消息
   - 详细的日志记录
   - 进度反馈

3. **错误处理**
   - 捕获具体异常
   - 提供有用信息
   - 优雅降级

### 设计模式

1. **工具函数模式**
   - 单一功能
   - 无状态
   - 易于测试

2. **策略模式**(改进版)
   - 支持多种转换策略
   - 可插拔的处理器

3. **建造者模式**(改进版)
   - 链式调用
   - 灵活配置

## 总结

这是一个**典型的Python工具函数**，体现了**简单、清晰、实用**的设计哲学。

**优点**:
- 代码简洁，易于理解
- 功能单一，职责明确
- 错误处理基本完整
- 用户体验友好

**不足**:
- 功能过于简单
- 缺乏配置灵活性
- 扩展性有限
- 大文件处理能力不足

**学习价值**:
- 掌握Python工具函数的编写规范
- 学习如何设计健壮的错误处理
- 理解类型注解和文档字符串的重要性
- 认识简单设计的价值

**改进方向**:
- 支持更多输入输出格式
- 添加配置和批处理能力
- 优化大文件处理
- 增强错误恢复

**核心启示**: 并不是所有代码都需要复杂的架构，有时一个清晰的函数就能解决大部分问题。这个脚本展示了**最小可行产品(MVP)**的设计思路——先实现核心功能，再逐步演进。

在生物信息学领域，这类简单的工具函数非常常见，它们往往被作为更大工作流的组成部分使用，因此**简洁性和可靠性**比**功能丰富性**更重要。
