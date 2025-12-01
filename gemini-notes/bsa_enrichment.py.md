# 代码结构优化-bsa_enrichment

本篇分析 `rnaseq/bsa_enrichment.py` 脚本，该脚本用于进行 GO 和 KEGG 富集分析。脚本虽然功能完整，但在路径管理、异常处理和可视化模块方面有待改进。

## 1. 路径管理：配置化与对象化
### 📌 问题（Line 157–162）
在 `bsa_enrichment` 函数中，硬编码了多个文件名拼接逻辑。如果文件命名规则发生变化，需要修改多处代码。

```python
    kegg_id_map = ann_dir / "kegg.idmap.csv"
    kegg_name = ann_dir / "kegg.name.csv"
    go_name = ann_dir / "go.name.csv"
    # ...
```

### ✅ 改进：使用配置对象或资源管理器
定义一个 `AnnotationResources` 类来管理所有相关文件的路径。

```python
@dataclass
class AnnotationResources:
    base_dir: Path
    
    @property
    def kegg_id_map(self) -> Path:
        return self.base_dir / "kegg.idmap.csv"
    
    @property
    def go_name(self) -> Path:
        return self.base_dir / "go.name.csv"
    
    # ... 其他属性
    
    def validate(self):
        # 检查所有文件是否存在
        if not self.kegg_id_map.exists():
            raise FileNotFoundError(f"Missing {self.kegg_id_map}")
```

## 2. 可视化模块：解耦绘图逻辑
### 📌 问题（Line 61–82）
`plot` 函数内部混合了数据预处理（计算 logPvalue）、样式设置（seaborn）和文件保存。且硬编码了图片尺寸计算逻辑，难以适应不同数量级的条目。

```python
    plot_height = 2 + len(plot_data) * 0.2
    text_length = plot_data.Description.map(lambda x: len(x)).max()
    plot_width = 6 + text_length * 0.1
```

### ✅ 改进：独立的绘图类
创建一个 `EnrichmentPlotter` 类，允许自定义绘图参数。

```python
class EnrichmentPlotter:
    def __init__(self, data: pd.DataFrame):
        self.data = data

    def preprocess(self, top_n=30):
        # 数据处理逻辑
        pass

    def plot_barplot(self, ax=None, **kwargs):
        # 纯粹的绘图逻辑，接受 matplotlib axes 对象
        pass

    def save(self, path: Path, width=None, height=None):
        # 自动计算或使用指定尺寸保存
        pass
```

## 3. 异常处理：装饰器模式
### 📌 问题（Line 174–196）
在主循环中，对 GO 和 KEGG 分析分别进行了重复的 `try-except` 块编写，代码冗余。

```python
            try:
                go(...)
            except Exception as e:
                logger.error(...)
            
            try:
                kegg(...)
            except Exception as e:
                logger.error(...)
```

### ✅ 改进：使用装饰器或上下文管理器
定义一个处理异常的装饰器，或者提取通用的执行函数。

```python
def safe_execute(func, description, *args, **kwargs):
    logger.info(f"Running {description}...")
    try:
        func(*args, **kwargs)
    except Exception as e:
        logger.error(f"Error running {description}")
        logger.error(e)

# 调用
safe_execute(go, "GO enrichment", gene_df=..., ...)
safe_execute(kegg, "KEGG enrichment", gene_df=..., ...)
```

## 🎯 改进总结
### 提升点
*   **健壮性**：资源管理类确保在开始分析前所有依赖文件都存在。
*   **灵活性**：绘图模块解耦后，可以轻松调整样式或集成到其他报告中。
*   **整洁性**：消除重复的异常处理代码，使主逻辑更清晰。

### 适用场景
*   批量处理生物信息学分析流程。
*   需要生成高质量出版级图片的脚本。
