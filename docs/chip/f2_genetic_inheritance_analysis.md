# F2代遗传一致性分析方法与参数设计

本设计文档描述了针对 `chip/genetic_inheritance.py` 工具从 F1 代遗传一致性检查扩展到支持 F2 代检查的方法学。**该功能已实现**。

## 1. 背景

现有的 `genetic_inheritance.py` 工具主要设计用于 F1 代的一致性验证（例如：杂交子代）。
在 F1 模型下，子代的每一个位点上的两个等位基因，必须严格地：
- 一个来自父本（Father）
- 一个来自母本（Mother）

但在实际的育种或群体遗传学研究中，经常会遇到 F2 代（由 F1 代个体间相互交配或自交产生的后代）。
在 F2 代中，子代的两个等位基因只要存在于原始双亲的"等位基因池"中即可。换句话说，由于在减数分裂和配子形成过程中的重组和自由组合，F2 代可能出现亲本所有的等位基因的任意两两组合，包括重新回到纯合状态（如具有和祖父本或祖母本完全一样的纯合基因型）。

## 2. 遗传模型设计

假设一个变异位点上：
- 亲本 1 的基因型等位基因为：$P1 = \{a, b\}$
- 亲本 2 的基因型等位基因为：$P2 = \{c, d\}$

### 2.1 F1 模型 (默认行为)
在 F1 代中，子代必然结合来自 P1 的一个等位基因和来自 P2 的一个等位基因。
因此 F1 代的有效基因型集合为（不考虑组合顺序）：
$$ F1_{expected} = \{(x, y) \mid x \in P1 \text{ 且 } y \in P2 \} $$
即：$\{(a,c), (a,d), (b,c), (b,d)\}$。
如果子代基因型不在此集合内，则视为遗传不一致。

### 2.2 F2 模型 (`--mode f2`)
在 F2 代中，它的亲本是 F1 个体。而所有 F1 个体等位基因的来源，追溯回去都是 P1 和 P2 提供的四个等位基因（即 $a, b, c, d$）。
因此，在理论上，F2 代个体可能获得的等位基因"池"（Pool）为亲本等位基因的并集（去重）：
$$ Pool = P1 \cup P2 = \{a, b, c, d\} $$

那么，F2 代个体的任何一个合法的基因型，必须由这个 Pool 中的任意两个等位基因组合而成。
F2 代的有效基因型集合为：
$$ F2_{expected} = \{(x, y) \mid x \in Pool \text{ 且 } y \in Pool \} $$

> **特例：**
> 对于普通的二倍体双亲纯合系杂交：
> 父本(AA) × 母本(BB) $\rightarrow$ F1(AB)
> - F1模型检查子代是否为 (A, B)
> - F2模型检查子代等位基因是否在集合 {A, B} 内，因此子代合法基因型为 (A,A), (A,B), (B,B)。

### 2.3 适用边界（重要）
当前实现的 `--mode f2` 是**宽松基因池模型**，用于遗传一致性质控放宽：
- 判定标准是“子代每个等位基因都在父母等位基因并集中”。
- 它不等价于“严格的F2家系重建模型”，不会基于具体F1来源推导完整重组概率与可达组合。

## 3. 使用方法

### 3.1 CLI 参数

```bash
# F1模式测试(默认)
python chip/genetic_inheritance.py test.vcf family.tsv output_f1.csv

# F2模式测试
python chip/genetic_inheritance.py test.vcf family.tsv output_f2.csv --mode f2

# 查看帮助
python chip/genetic_inheritance.py --help
```

### 3.2 参数说明

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--mode` | `f1` | 遗传一致性检查模式: `f1`=直接亲子(F1代), `f2`=原始亲本到F2代 |

### 3.3 行为对比

| 场景 | F1验证 | F2验证 |
|------|--------|--------|
| 父: A/A, 母: B/B, 子: A/B | 通过 | 通过 |
| 父: A/A, 母: B/B, 子: A/A | 不通过 | 通过 |
| 父: A/B, 母: C/D, 子: A/C | 通过 | 通过 |
| 父: A/B, 母: C/D, 子: A/E | 不通过 | 不通过 |

## 4. 代码实现

### 4.1 CLI 参数定义
```python
def main(
    vcf_file: Path,
    family_file: Path,
    output_file: Path,
    max_variants: Optional[int] = None,
    max_alleles: int = 2,
    mode: str = typer.Option(
        "f1",
        "--mode",
        help="遗传一致性检查模式: f1=直接亲子(F1代), f2=原始亲本到F2代",
    ),
):
    """主程序"""
    analyzer = VCFGeneticAnalyzer(mode=mode)
```

### 4.2 验证逻辑实现
```python
def check_variant_consistency(
    self, variant_id: str, father_gt: str, mother_gt: str, child_gt: str
) -> Tuple[bool, str]:
    """检查单个位点遗传一致性，支持F1和F2两种模式"""
    # ... 前置检查代码（解析基因型、检查缺失、二倍体检查）...

    if self.mode == "f1":
        # F1模式: 子代必须从父、母各继承一个等位基因
        possible_children = set()
        for f in father_alleles:
            for m in mother_alleles:
                possible_children.add(tuple(sorted([f, m])))

        child_tuple = tuple(sorted(child_alleles))

        if child_tuple in possible_children:
            return True, "遗传一致"
        else:
            return (
                False,
                f"子代({child_alleles})不可能由父({father_alleles})母({mother_alleles})组合得到",
            )
    else:
        # F2模式: 子代两个等位基因只需存在于基因池(父∪母)中
        gene_pool = set(father_alleles) | set(mother_alleles)

        child_in_pool = [allele in gene_pool for allele in child_alleles]

        if all(child_in_pool):
            return True, "遗传一致"
        else:
            missing_alleles = [a for a, in_pool in zip(child_alleles, child_in_pool) if not in_pool]
            return (
                False,
                f"子代等位基因{missing_alleles}不在基因池({list(gene_pool)})中"
            )
```

## 5. 向后兼容性

1. **默认行为不变**：`--mode` 默认值为 `f1`，现有命令无需修改即可正常运行
2. **输出格式不变**：所有输出文件格式与之前完全一致
3. **错误信息保持中文**：新增的错误信息同样使用中文

## 6. 总结

通过 `--mode` 参数，工具现已支持 F1 严格模式与 F2 宽松基因池模式的遗传一致性质控。
