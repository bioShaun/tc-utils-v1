# 代码结构优化-finger_snp

本篇分析 `chip/finger_snp.py` 脚本，该脚本用于选择最小 SNP 集合以区分个体。脚本实现了多种算法（贪心、暴力、遗传），但设计模式上可以进一步解耦，以提高可扩展性和可测试性。

## 1. 算法解耦：策略模式（Strategy Pattern）
### 📌 问题（Line 92–318）
`SNPMinimizer` 类中直接实现了三种不同的算法 (`greedy_selection`, `brute_force_search`, `genetic_algorithm`)。如果未来要添加新算法（如模拟退火），需要修改这个类，违反了开闭原则（Open-Closed Principle）。

```python
class SNPMinimizer:
    def greedy_selection(self): ...
    def brute_force_search(self): ...
    def genetic_algorithm(self): ...
```

### ✅ 改进：定义算法接口
使用策略模式将每种算法封装为独立的类。

```python
from abc import ABC, abstractmethod
from typing import List, Tuple

class SelectionStrategy(ABC):
    @abstractmethod
    def select(self, data: np.ndarray) -> List[int]:
        pass

class GreedyStrategy(SelectionStrategy):
    def select(self, data: np.ndarray) -> List[int]:
        # 实现贪心算法
        pass

class GeneticStrategy(SelectionStrategy):
    def select(self, data: np.ndarray) -> List[int]:
        # 实现遗传算法
        pass

class SNPMinimizer:
    def __init__(self, strategy: SelectionStrategy):
        self.strategy = strategy

    def run(self, data):
        return self.strategy.select(data)
```

## 2. 职责分离：计算与展示解耦
### 📌 问题（Line 319–345）
`analyze_results` 方法混合了结果验证、数据提取和控制台打印。这使得该逻辑难以在 Web 服务或 GUI 中复用。

```python
    def analyze_results(self, selected_snps, selected_names):
        if not selected_snps:
            print("未找到有效解决方案") # 硬编码打印
            return
        # ... 打印更多信息
```

### ✅ 改进：返回结构化结果
计算类只负责返回数据对象，展示逻辑由调用方处理。

```python
@dataclass
class SelectionResult:
    selected_indices: List[int]
    selected_ids: List[str]
    is_valid: bool
    execution_time: float
    # ... 其他元数据

class SNPMinimizer:
    def run(self) -> SelectionResult:
        # ... 计算逻辑
        return SelectionResult(...)

# 调用方负责展示
def print_report(result: SelectionResult):
    print(f"选择的SNP数量: {len(result.selected_indices)}")
```

## 3. 性能优化：向量化计算
### 📌 问题（Line 131–144）
`_count_distinguished_pairs` 使用双重循环遍历所有个体对，时间复杂度为 O(N^2)。对于大规模群体，这将非常慢。

```python
        for i in range(self.n_individuals):
            for j in range(i + 1, self.n_individuals):
                if not np.array_equal(selected_data[i], selected_data[j]):
                    distinguished_pairs += 1
```

### ✅ 改进：利用 NumPy 广播或哈希
使用 NumPy 的广播机制或将行转换为哈希值进行快速比较。

```python
    def _can_distinguish_all(self, snp_indices: List[int]) -> bool:
        selected_data = self.snp_data[:, snp_indices]
        # 利用 unique 快速检查唯一行数
        # axis=0 表示按行比较
        unique_rows = np.unique(selected_data, axis=0)
        return len(unique_rows) == self.n_individuals
```

## 🎯 改进总结
### 提升点
*   **可扩展性**：新增算法只需添加一个新的策略类，无需修改现有代码。
*   **复用性**：核心逻辑不再依赖 `print`，可轻松集成到 Web API 或其他工具链中。
*   **性能**：利用 NumPy 特性替代 Python 循环，显著提升在大数据集上的运行速度。

### 适用场景
*   算法比较与基准测试。
*   需要高性能计算的数据分析工具。
