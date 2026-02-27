# F2代遗传一致性分析方法与参数设计

本设计文档描述了针对 `chip/genetic_inheritance.py` 工具从 F1 代遗传一致性检查扩展到支持 F2 代检查的方法学。**用户要求目前只需输出分析方法文档，不需要修改代码**。

## 1. 背景

现有的 `genetic_inheritance.py` 工具主要设计用于 F1 代的一致性验证（例如：杂交子代）。
在 F1 模型下，子代的每一个位点上的两个等位基因，必须严格地：
- 一个来自父本（Father）
- 一个来自母本（Mother）

但在实际的育种或群体遗传学研究中，经常会遇到 F2 代（由 F1 代个体间相互交配或自交产生的后代）。
在 F2 代中，子代的两个等位基因只要存在于原始双亲的“等位基因池”中即可。换句话说，由于在减数分裂和配子形成过程中的重组和自由组合，F2 代可能出现亲本所有的等位基因的任意两两组合，包括重新回到纯合状态（如具有和祖父本或祖母本完全一样的纯合基因型）。

## 2. 遗传模型设计

假设一个变异位点上：
- 亲本 1 的基因型等位基因为：$P1 = \{a, b\}$
- 亲本 2 的基因型等位基因为：$P2 = \{c, d\}$

### 2.1 现有的 F1 模型 (默认行为)
在 F1 代中，子代必然结合来自 P1 的一个等位基因和来自 P2 的一个等位基因。
因此 F1 代的有效基因型集合为（不考虑组合顺序）：
$$ F1_{expected} = \{(x, y) \mid x \in P1 \text{ 且 } y \in P2 \} $$
即：$\{(a,c), (a,d), (b,c), (b,d)\}$。
如果子代基因型不在此集合内，则视为遗传不一致。

### 2.2 扩展的 F2 模型 (`--pop f2`)
在 F2 代中，它的亲本是 F1 个体。而所有 F1 个体等位基因的来源，追溯回去都是 P1 和 P2 提供的四个等位基因（即 $a, b, c, d$）。
因此，在理论上，F2 代个体可能获得的等位基因“池”（Pool）为亲本等位基因的并集（去重）：
$$ Pool = P1 \cup P2 = \{a, b, c, d\} $$

那么，F2 代个体的任何一个合法的基因型，必须由这个 Pool 中的任意两个等位基因组合而成。
F2 代的有效基因型集合为：
$$ F2_{expected} = \{(x, y) \mid x \in Pool \text{ 且 } y \in Pool \} $$

> **特例：**
> 对于普通的二倍体双亲纯合系杂交：
> 父本(AA) × 母本(BB) $\rightarrow$ F1(AB)
> - F1模型检查子代是否为 (A, B)
> - F2模型检查子代等位基因是否在集合 {A, B} 内，因此子代合法基因型为 (A,A), (A,B), (B,B)。

## 3. 代码实现指导（未来实施参考）

如果未来需要在 `genetic_inheritance.py` 中实现该功能，可以参考以下设计：

### 3.1 CLI 参数调整
在 `main()` 函数或 CLI 入口添加群体类型参数：
```python
import typer
from typing_extensions import Annotated

def main(
    # ... 其他参数
    pop: Annotated[str, typer.Option("--pop", help="群体类型，支持f1(默认)或f2")] = "f1"
):
    # 将 pop 传入 analyzer 对象并保存为类属性，供后续验证时使用
    analyzer.pop = pop.lower()
```

### 3.2 验证逻辑重构
在 `check_variant_consistency()` 方法中，针对 `pop` 类型采取不同的校验策略：

```python
def check_variant_consistency(self, variant_id: str, father_gt: str, mother_gt: str, child_gt: str) -> Tuple[bool, str]:
    # 1. 获取并解析真实的等位基因
    father_alleles = list(self.resolve_alleles(variant_id, father_indices))
    mother_alleles = list(self.resolve_alleles(variant_id, mother_indices))
    child_alleles = list(self.resolve_alleles(variant_id, child_indices))
    
    # ... (省略缺失检查和二倍体检查等预处理) ...

    possible_children = set()
    
    if getattr(self, "pop", "f1") == "f2":
        # === F2代验证逻辑 ===
        # 等位基因池是父母双方提供的所有等位基因的并集
        parent_allele_pool = list(set(father_alleles + mother_alleles))
        
        # 子代的两个等位基因可以从 pool 中任意取两个组成
        for a1 in parent_allele_pool:
            for a2 in parent_allele_pool:
                possible_children.add(tuple(sorted([a1, a2])))
                
        child_tuple = tuple(sorted(child_alleles))
        
        if child_tuple in possible_children:
            return True, "遗传一致"
        else:
            return False, f"F2子代({child_alleles})的等位基因不在祖代池({parent_allele_pool})中"
            
    else:
        # === 原有的 F1 代验证逻辑 ===
        for f in father_alleles:
            for m in mother_alleles:
                possible_children.add(tuple(sorted([f, m])))
                
        child_tuple = tuple(sorted(child_alleles))

        if child_tuple in possible_children:
            return True, "遗传一致"
        else:
            return False, f"F1子代({child_alleles})不可能由父({father_alleles})母({mother_alleles})组合得到"
```

## 4. 总结
通过扩展参数和基于群体类型的遗传重组算法判断，即可灵活地支持 F1/F2 及不同群体世代的遗传一致性质控。当前暂不修改代码，待需求明确实施时，可按此文档指引直接接入。
