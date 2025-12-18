# `replace-probe-by-distance-codex.py`

This script replaces a given set of probes with optimal candidates from a larger pool based on genomic distance and other metrics. It is designed to find the best possible replacement for a probe within a specified maximum distance.

## 核心替换逻辑

脚本为每个待替换的探针（probe）寻找最佳的候选探针。搜索策略如下：

1.  **最大范围**: 只在指定的最大距离内（`--max-distance-kb`）搜索同一染色体上的候选探针。
2.  **优先级（Priority）**: 优先选择 `priority` 值最低的候选者。
3.  **距离窗口**: 在同一 `priority` 级别下，根据距离远近将候选者分组到不同的窗口（`--window-kb`）。优先选择距离原始位置更近的窗口。
4.  **MAF**: 在同一窗口内，选择 `maf`（Minor Allele Frequency）最高的候选者。
5.  **精确距离**: 如果以上条件都相同，则选择实际距离（`distance_bp`）最近的。
6.  **唯一性**: 每个候选探针只能被使用一次。

## 使用方法

```bash
python replace-probe-by-distance-codex.py [OPTIONS] REPLACE_FILE CANDIDATE_FILE OUT_FILE
```

## 参数说明

### Positional Arguments

-   `REPLACE_FILE`: (必须) 一个表格文件，列出了需要被替换的探针。
-   `CANDIDATE_FILE`: (必须) 一个表格文件，包含了所有可用于替换的候选探针。
-   `OUT_FILE`: (必须) 输出文件的路径，用于保存替换后的探针表格。可以是 `.tsv` 或 `.xlsx` 格式。

### Options

-   `--window-kb, -w INTEGER`: 搜索窗口的大小（单位：kb）。默认为 `10`。
-   `--max-distance-kb, -m INTEGER`: 允许搜索的最大距离（单位：kb）。默认为 `100`。
-   `--mapping, -p FILE`: 指定一个文件路径，用于输出原始探针与替换探针之间的映射关系。如果未提供，则默认在 `OUT_FILE` 同目录下生成一个 `_mapping.tsv` 文件。
-   `--allow-missing`: 一个开关选项。如果设置了此项，即使某些探针在最大范围内找不到任何可用的替换，脚本也会继续运行并生成部分结果。同时，会额外输出两个文件：
    -   `*_no_replacement.tsv`: 无法找到替换的探针列表。
    -   `*_remaining_candidates.tsv`: 未被使用的剩余候选探针列表。
-   `--help`: 显示帮助信息。

## 输入文件格式

输入文件应为制表符分隔的文本文件（`.tsv`, `.txt`）或类似的表格格式。

### `REPLACE_FILE` 必需列

-   `chrom`: 染色体名称
-   `pos`: 探针位置
-   `id`: 探针ID
-   `target_id`: 目标ID

### `CANDIDATE_FILE` 必需列

-   `chrom`: 染色体名称
-   `pos`: 探针位置
-   `id`: 探针ID
-   `target_id`: 目标ID
-   `maf`: Minor Allele Frequency (等位基因频率)
-   `priority`: 优先级 (数值越小，优先级越高)

## 输出文件

-   **`OUT_FILE`**: 包含最终被选中的替换探针。其列信息来自候选探针文件，并额外增加了以下几列来说明其来源：
    -   `origin_chrom`: 原始探针的染色体
    -   `origin_pos`: 原始探针的位置
    -   `distance_to_origin_bp`: 与原始探针的距离（碱基对）

-   **`_mapping.tsv`**: 记录了原始探针和替换探针的详细对应关系，方便追溯。

## 使用示例

```bash
# 基本用法
python draft/replace-probe-by-distance-codex.py \
    probes_to_replace.tsv \
    candidate_probes.tsv \
    replaced_probes.tsv

# 指定搜索范围并允许部分失败
python draft/replace-probe-by-distance-codex.py \
    --window-kb 5 \
    --max-distance-kb 50 \
    --allow-missing \
    --mapping replacement_map.tsv \
    probes_to_replace.tsv \
    candidate_probes.tsv \
    replaced_probes_partial.xlsx
```
