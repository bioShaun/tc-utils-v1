# cdsCovEvaluation_bamdst_02x_polars.py - Polars数据分析实践

## 代码概览

**文件路径**: panel/cdsCovEvaluation_bamdst_02x_polars.py

**功能描述**: 使用Polars库分析多个样本的CDS区域覆盖度，计算每个区域的0.2x覆盖比例。

**核心技术**: Polars DataFrame操作、自定义异常体系、loguru日志系统、输入验证。

## 代码结构分析

### 架构设计

```
输入验证层 (validate_input_path, validate_dataframe)
        ↓
数据加载层 (load_bed_files, load_single_depth_file)
        ↓
数据转换层 (merge_chr, calculate_coverage_ratio)
        ↓
输出层 (write_output)
```

### 优点

1. **完整的异常体系** (第25-35行)
   ```python
   class CoverageAnalysisError(Exception):
       pass

   class DataValidationError(CoverageAnalysisError):
       pass

   class FileProcessingError(CoverageAnalysisError):
       pass
   ```
   - 层次化异常设计，便于调用方处理
   - 自定义异常类型清晰表达错误语义
   - 使用继承构建错误分类体系

2. **健壮的输入验证**
   ```python
   def validate_input_path(path: Path, path_type: str = "directory") -> None:
       if not path.exists():
           raise DataValidationError(f"{path_type.capitalize()} does not exist: {path}")
   ```
   - 统一验证函数，避免重复代码
   - 早期失败原则(fail-fast)
   - 清晰的错误消息

3. **Polars高效操作**
   - **水平拼接** (第194行): `pl.concat(df_list, how="horizontal")`
   - **链式调用** (第156-164行): 使用`.`操作符链式调用
   - **向量化计算** (第212-216行): `sum_horizontal()`一次计算所有列的和
   - **惰性求值**(部分使用): 减少中间结果创建

4. **日志系统设计** (第12-19行)
   ```python
   class InterceptHandler(logging.Handler):
       def emit(self, record):
           logging.getLogger(record.name).handle(record)
   ```
   - 桥接loguru和标准logging
   - 统一日志格式和输出

### 关键技术点

1. **坐标变换** (第60-97行)
   ```python
   merged_df = df.join(split_bed_df, on="chrom", how="inner")
   merged_df = merged_df.with_columns(
       (pl.col("start") + pl.col("offset")).alias("start"),
       (pl.col("end") + pl.col("offset")).alias("end"),
   )
   ```
   - 使用DataFrame连接实现染色体坐标转换
   - 通过`with_columns`同时计算多个新列

2. **动态列选择** (第82-86行)
   ```python
   final_cols = ["new_chrom", "start", "end"] + [
       col for col in merged_df.columns
       if col not in ["new_chrom", "start", "end", "offset", "offset_end", "chrom"]
   ]
   ```
   - 动态构建列列表，支持任意数量的附加列
   - 避免硬编码列名

3. **批处理模式** (第176-203行)
   - 遍历目录中的所有depth文件
   - 按需筛选样本列表
   - 构建BED矩阵进行批量分析

### 可维护性分析

#### 优秀实践

1. **类型注解完整**
   - 所有函数参数和返回值都有类型提示
   - 使用`Optional`、`Union`处理可选参数

2. **函数职责单一**
   - 每个函数只负责一个明确的任务
   - 易于测试和复用

3. **配置参数集中**
   - `DEFAULT_DEPTH_THRESHOLD`、`DEFAULT_FLOAT_PRECISION`等常量
   - 便于调整算法参数

#### 潜在问题

1. **错误恢复机制缺失**
   ```python
   # 现状：单个文件失败则整个分析失败
   for bed_file in bed_list:
       df_i = load_single_depth_file(bed_file, sample_name, chrom_prefix)
       if df_i is not None:
           df_list.append(df_i)
   ```
   - **影响**: 部分样本损坏导致整个分析中断
   - **建议**: 支持跳过错误样本并继续处理

2. **配置硬编码**
   ```python
   # 深度阈值硬编码
   depth_threshold = mean_depth * DEFAULT_DEPTH_THRESHOLD
   ```
   - **影响**: 难以适应不同项目需求
   - **建议**: 支持命令行参数配置

3. **日志级别混淆**
   ```python
   # 使用标准logging但输出格式与loguru不一致
   logging.info("Starting CDS coverage analysis")
   logger.info("Processing sample: {sample_name}", sample_name=sample_name)
   ```
   - **影响**: 日志输出格式不统一
   - **建议**: 统一使用loguru或标准logging

## Polars最佳实践

### 1. 数据加载

```python
# 多文件批量读取
df_list = []
for bed_file in bed_list:
    df_i = load_single_depth_file(bed_file, sample_name)
    if df_i is not None:
        df_list.append(df_i)

# 水平拼接构建样本矩阵
df_matrix = pl.concat(df_list, how="horizontal")
```

### 2. 数据清洗

```python
# 过滤空值和无效数据
df = df.filter(pl.col("扣费时间").is_not_null())

# 类型转换与填充
df = df.with_columns(
    pl.col("CPU核*时").cast(pl.Float64, strict=False).fill_null(0),
)
```

### 3. 聚合计算

```python
# 分组聚合
monthly_stats = df.group_by(["年月", "退出状态"]).agg([
    pl.sum("CPU核*时").alias("CPU核*时"),
    pl.sum("金额(元)").alias("金额(元)"),
])
```

### 4. 列操作

```python
# 计算覆盖比例
cover_sum = df_matrix.sum_horizontal()
cover_ratio = cover_sum / df_matrix.width

# 计算新坐标
merged_df = merged_df.with_columns(
    (pl.col("start") + pl.col("offset")).alias("start"),
    (pl.col("end") + pl.col("offset")).alias("end"),
)
```

## 性能特点

### Polars vs Pandas

1. **惰性求值**
   - Polars使用查询优化器，最小化内存使用
   - Pandas总是立即执行计算

2. **并行执行**
   - Polars自动并行化操作
   - Pandas需要手动指定

3. **内存效率**
   - Polars使用Apache Arrow后端，内存效率更高
   - 更少的数据拷贝

### 性能优化建议

1. **避免`to_pandas()`转换**
   - 现状: 第171行使用`to_pandas()`后写入Excel
   - 建议: 使用Polars原生`write_excel()`(如果支持)或`write_ipc()`

2. **使用IPC格式缓存**
   - 参考`hpc_cost_calculator_polars.py`的做法
   - 避免重复读取大文件

3. **惰性模式**
   ```python
   # 使用lazy API处理大数据集
   df = pl.scan_csv("large_file.csv")
   result = df.filter(...).group_by(...).agg(...).collect()
   ```

## 学习要点

### Python技巧

1. **异常层次设计**
   - 顶层通用异常 → 细分异常类型
   - 便于调用方精确捕获和处理

2. **类型系统**
   - `Optional[Type]`: 可能为None的类型
   - `Union[Type1, Type2]`: 多种可能类型
   - 提升代码可读性和IDE支持

3. **日志系统整合**
   - `InterceptHandler`桥接不同日志库
   - 统一日志格式和输出目标

### 数据处理模式

1. **验证-转换-输出**
   - 先验证输入有效性
   - 转换数据格式和内容
   - 输出到目标文件

2. **批处理模式**
   - 遍历多个输入文件
   - 批量处理和合并
   - 生成统一输出

## 改进建议

### 1. 配置管理

```python
@dataclass
class CoverageConfig:
    """覆盖度分析配置"""
    depth_threshold: float = DEFAULT_DEPTH_THRESHOLD
    float_precision: int = DEFAULT_FLOAT_PRECISION
    chrom_prefix: Optional[str] = None
    log_level: str = "INFO"

    @staticmethod
    def from_typer(ctx: typer.Context):
        """从Typer上下文创建配置"""
        return CoverageConfig(
            depth_threshold=ctx.params.get('depth_threshold', DEFAULT_DEPTH_THRESHOLD),
            # ...
        )
```

### 2. 错误恢复

```python
def load_bed_files_robust(bed_dir: Path) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """支持跳过错误样本的版本"""
    bed_df, df_matrix = load_bed_files(bed_dir)
    if df_matrix.is_empty():
        logger.warning("部分样本加载失败，返回可用数据")
    return bed_df, df_matrix
```

### 3. 进度追踪

```python
from tqdm import tqdm

for bed_file in tqdm(bed_list, desc="加载深度文件"):
    df_i = load_single_depth_file(bed_file, sample_name)
    # ...
```

## 测试策略

1. **单元测试**
   - 输入验证函数
   - 数据转换函数
   - 计算函数

2. **集成测试**
   - 完整流程测试
   - 错误处理测试

3. **性能测试**
   - 大文件处理能力
   - 内存使用监控

## 总结

这是一个**现代Python数据处理的优秀范例**，展示了：

1. **Polars的高效DataFrame操作**
2. **完整的异常处理体系**
3. **健壮的输入验证机制**
4. **清晰的功能分离**

**值得学习的地方**:
- 异常层次的设计思维
- Polars链式调用的优雅写法
- 类型注解的全面使用
- 输入验证的防御性编程

**改进方向**:
- 统一日志系统
- 增加错误恢复机制
- 配置外部化
- 添加进度追踪

整体而言，这是一个**生产就绪的高质量脚本**，体现了作者对Python数据处理和工程实践的深刻理解。
