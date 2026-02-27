# AGENTS.md - tc-pytools 代理协作约束（v4）

版本日期：2026-02-27

## 目的（只写防错信息）

本文件只保留“代理容易做错的点”和“必须遵守的护栏”。
不要在这里维护目录树、脚本清单或固定构建/测试命令；这些信息应从仓库现状（`pyproject.toml`、README、代码）自动发现。

## 必须遵守（MUST）

- 先探索再决策：修改或新增前，先搜索并阅读真实入口与现有实现模式，禁止凭猜测创建路径/模块。
- CLI 与核心逻辑分离：业务逻辑必须可测试；CLI 只负责参数、I/O、退出码。
- 日志与输出：禁止散落调试 `print()`；优先 `loguru`；CLI 最终用户输出可用 `typer.echo()` 或 `rich`。
- 错误处理：禁止裸 `except:`；失败时必须 `raise` 或 `typer.Exit(code=1)`，并包含上下文（输入路径、关键参数、记录定位信息）。
- 大文件处理：VCF/FASTA/BAM/超大表格禁止一次性读入内存，优先流式、分块或惰性计算。
- 路径处理：使用 `pathlib.Path`；禁止硬编码集群绝对路径（如 `/mnt/...`）。

## CLI 护栏（MUST + 短模板）

- 新 CLI 优先 Typer + `Annotated`。
- `--help` 必须包含：功能简述、使用示例、输出格式；参数必须有中文 help。
- 使用短模板骨架，不在本文件维护长脚手架。

```python
from inspect import cleandoc
from pathlib import Path
from typing import Annotated

import typer

MODULE_HELP = cleandoc(
    """
    脚本功能简述。

    \b
    使用示例:
      python script.py in.tsv out.tsv

    \b
    输出格式:
      - out.tsv: 列说明
    """
)

app = typer.Typer(help=MODULE_HELP, no_args_is_help=True)


@app.command(help=MODULE_HELP)
def main(
    input_file: Annotated[Path, typer.Argument(help="输入文件")],
    output_file: Annotated[Path, typer.Argument(help="输出文件")],
) -> None:
    ...
```

## 数据处理建议（SHOULD）

- 中大规模表格 SHOULD 优先 Polars；小规模或兼容场景可用 pandas。
- 若尝试 Polars 但环境不满足，必须显式报错或明确降级策略，禁止静默改变结果语义。
- 生信格式优先使用专业库：VCF 用 `cyvcf2`/`pysam`，FASTA 用 `pyfaidx`，BAM/SAM 用 `pysam`。

## Gotchas（强制保留）

- 仓库历史包袱较重，存在多种代码风格并存；旧风格脚本不应直接作为新功能模板。
- 遗留风格示例：`bsa/bsa_region.py`、`fq/check_q30.py`（argparse/pandas 传统写法）。
- 新功能参考示例：`gtf/split_large_genome_bed_v2.py`、`chip/panel2bed.py`、`exome/splitBed_v2.py`（更现代的 CLI 与工程化写法）。
- 对 `draft/` 下脚本默认视作实验性实现，复用前必须先核验边界条件与测试覆盖。

## 遇到惊吓请立刻报告（STOP 条款）

如果遇到任何让你惊讶或困惑的地方（例如文档与代码不一致、入口不明确、依赖冲突、目录组织反直觉），请立刻：

1. 停止继续“猜测式修改”（不要凭空新建文件来凑通）
2. 明确指出困惑点（文件路径、期望行为、实际观察）
3. 给出 1-2 个可选修复方向（改代码、改文档、补配置）并说明推荐方案
