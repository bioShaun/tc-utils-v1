# 代码结构优化-gwas report

本篇分析 `gwas/gwas-report.py` 脚本，该脚本使用 `reportlab` 生成 PDF 报告。虽然功能实现完整，但在代码结构、样式管理和可扩展性方面有较大的优化空间。

## 1. 样式管理：解耦样式与业务逻辑
### 📌 问题（Line 34–78）
在 `Graphs` 类的静态方法中，每次绘制元素时都重新获取样式表并修改属性。这种方式导致样式定义分散，难以统一管理和复用。

```python
    # 绘制标题
    @staticmethod
    def draw_title(title: str):
        # 获取所有样式表
        style = getSampleStyleSheet()
        # 拿到标题样式
        ct = style["Heading1"]
        # 单独设置样式相关属性
        ct.fontName = "NotoSansSC"  # 字体名
        ct.fontSize = 18  # 字体大小
        ...
```

### ✅ 改进：集中式样式管理器
创建一个单例或模块级的样式管理器，统一配置所有样式。

```python
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle

class StyleManager:
    _instance = None
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance.styles = getSampleStyleSheet()
            cls._instance._init_custom_styles()
        return cls._instance

    def _init_custom_styles(self):
        self.styles.add(ParagraphStyle(
            name='CustomTitle',
            parent=self.styles['Heading1'],
            fontName='NotoSansSC',
            fontSize=18,
            leading=50,
            alignment=1,
            spaceAfter=20
        ))
        # ... 其他样式定义

    def get(self, name):
        return self.styles[name]
```

## 2. 数据结构：使用 DataClass 封装上下文
### 📌 问题（Line 147–156）
`report` 函数接收大量参数，导致函数签名冗长，且参数之间缺乏逻辑关联。

```python
def report(
    species: str,
    genome: str,
    probe_tag: str,
    data_path: Path,
    report_path: Path,
    maf: Optional[float] = 0.05,
    # ... 更多参数
) -> None:
```

### ✅ 改进：定义报告上下文类
使用 `dataclass` 封装报告所需的数据和配置。

```python
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class ReportConfig:
    maf: float = 0.05
    gc_low: float = 0.3
    gc_high: float = 0.6
    match_criteria: str = "30,1"

@dataclass
class ReportContext:
    species: str
    genome_version: str
    probe_count: str
    data_dir: Path
    output_path: Path
    config: ReportConfig
```

## 3. 逻辑解耦：构建者模式（Builder Pattern）
### 📌 问题（Line 147–219）
`report` 函数内部混合了数据处理、页面布局和内容生成。随着报告内容增加，该函数会变得极难维护。

### ✅ 改进：Report Builder 类
将报告的构建过程拆分为独立的步骤。

```python
class PDFReportBuilder:
    def __init__(self, context: ReportContext):
        self.ctx = context
        self.elements = []
        self.styles = StyleManager()

    def add_header(self):
        # 添加页眉逻辑
        pass

    def add_basic_info(self):
        title = Paragraph("芯片设计报告", self.styles.get('CustomTitle'))
        self.elements.append(title)
        # ... 添加表格逻辑

    def add_plots(self):
        # ... 添加图片逻辑
        pass

    def build(self):
        doc = SimpleDocTemplate(str(self.ctx.output_path))
        doc.build(self.elements)
```

## 🎯 改进总结
### 提升点
*   **可维护性**：样式统一管理，修改字体或间距只需改动一处。
*   **可读性**：`report` 函数不再是几百行的面条代码，而是清晰的步骤调用。
*   **可扩展性**：新增章节只需在 Builder 中添加一个方法，不影响其他部分。

### 适用场景
*   生成复杂的 PDF/Excel/HTML 报告。
*   需要统一 UI/VI 风格的项目。
