# TC PyTools

集合了一些在生信流程中常用的 Python/脚本工具，按业务场景分模块（panel、primer、fasta、vcf 等）。本仓库没有统一的入口，每个脚本独立使用，建议在对应目录下查阅专属说明。

## 使用方式
- 依赖安装：`pip install -r requirements.txt`
- 本地检查：`pre-commit install` 后，提交前会自动运行格式化（ruff-format）、静态检查（ruff、pyright）和 `pytest`
- 手动运行检查：`pre-commit run --all-files`
- 单元测试：`pytest`

## 文档链接
- primer 模块：`primer/README_add_flank_sequence.md`（介绍加引物侧翼序列脚本的用法）
- primer 模块：`primer/README.md`（align-kasp 等脚本说明）

提交新脚本时，建议在对应目录添加 README 描述输入/输出示例和依赖，并在这里补充链接，便于快速查阅。***
