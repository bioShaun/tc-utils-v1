# primer

本目录包含与引物/探针相关的脚本与测试。

## align-kasp.py

根据 KASP 引物表（FAM/VIC/Common 三条引物序列）调用 `blastn` 将引物比对到参考基因组（BLAST DB），并输出每条引物的基因组位置与 SNP 位点（3' 端锚定）。
为保证 SNP 位点的准确性，脚本要求比对结果覆盖引物的 3' 端（`qend` 接近引物长度）；若出现 3' 端轻微 clip（例如 `qend = length - 1`），会按链方向外推到引物真实 3' 端坐标。

### 依赖

- Python 依赖：见仓库根目录 `requirements.txt`（脚本使用 `pandas/numpy/typer/loguru` 等；`delegator` 为可选依赖，没有时会回退到 `subprocess`）
- 外部程序：`blastn`（NCBI BLAST+），以及已构建好的 BLAST 数据库（`makeblastdb`）

### 输入

KASP 表：无表头、4 列、TAB 分隔：

1. `name`：引物/位点名称（唯一标识）
2. `fam`：FAM 引物序列
3. `vic`：VIC 引物序列
4. `common`：Common 引物序列

序列允许 IUPAC 码（脚本会将每个 IUPAC 碱基映射为一个 ATGC 碱基用于比对）。

### 输出

在 `out_dir` 下生成：

- `kasp.fam.fa` / `kasp.vic.fa` / `kasp.common.fa`：用于比对的 FASTA
- 对应的 `*.blasttab.tsv`：`blastn -outfmt 6` 的结果（脚本会自动生成文件名）

在输入文件同目录下生成：

- `<input>.pos.tsv`：汇总表（包含 `chrom`、`snp_pos`、`primer_span` 等字段）

### 用法

```bash
python primer/align-kasp.py <kasp.tsv> <blast_db_prefix> <out_dir> --threads 16 --max-mismatch 3
```

示例：

```bash
python primer/align-kasp.py kasp.tsv genome_db out.align-kasp --threads 8 --force --max-mismatch 3
```

### 测试

仓库使用 `pytest`，本脚本测试在 `primer/tests/test_align_kasp.py`：

```bash
pytest -q primer/tests/test_align_kasp.py
```
