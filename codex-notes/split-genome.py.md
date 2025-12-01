# split-genome.py

## 功能与输入输出
- 目的：在基因间隙处拆分超长染色体，生成新的 split FASTA、对应的 split GFF，以及记录拆分区间的 BED。
- 输入：`genome_fa`、`genome_fai`、`gff`；输出：`split_cat_bed`(拆分/未拆分的染色体区间)、`gff.split.gff`、`genome_fa.split.fa`。

## 核心流程
- 读 FAI -> `compute_split_windows` 计算每条染色体可拆分窗口（按比例截取两端，大小由 `split_size` 决定）。
- 读 GFF -> `parse_gff_with_pandas` 自动选择 feature 类型（优先 gene/mRNA/...）或使用指定 feature；`calculate_gaps` 找到相邻 feature 的间隔。
- `select_split_candidates` 在窗口内挑最大 gap，按染色体一对一选择并校验是否满足 `min_gene_gap`。
- `generate_split_chr_bed` 生成拆分点的 BED 定义；`finalize_unsplit_chromosomes` 补上无需拆分的染色体。
- `generate_split_gff` 依照 BED 重映射坐标，输出拆分后的 GFF；`generate_split_genome` 用 pyfaidx 按 BED 输出序列。

## 观察与隐患
- `generate_split_genome` 的换行长度参数未用到，循环里写死 60；`line_length` 形参被忽略。
- 全量 `pd.read_csv` 读取 GFF，超大文件可能内存高；且不支持 .gz。
- `compute_split_windows` 的分割窗口计算比较魔法（0.3/比例截取），缺少注释或引用来源，难以调整。
- 拆分时 `Fasta` 按染色体随机访问，如果 BED 被错误生成会直接 `ValueError`；缺少对 BED 数据自身的校验（是否覆盖范围、是否有交叉）。

## 改进建议
- 修复 `generate_split_genome` 使用 `line_length` 形参，或删掉形参/常量，避免误导。
- 为关键步骤写文档/注释：`split_size`、窗口公式、默认 feature 选择策略。
- 增加输入校验：BED 是否与 FAI 染色体一致、是否重叠；GFF 是否按染色体排序（可在 `calculate_gaps` 前排序）。
- 支持流式/分块读取 GFF，或允许传入预过滤的 GFF；支持 gzip 输入。
- CLI 体验：为 Typer 参数添加 help/类型限制，提供 `--feature-type`、`--line-length` 等选项的说明。

## 测试思路
- 构造小型 FASTA/FAI/GFF：
  - 单染色体、可拆分且 gap 充足 -> 生成 a/b、GFF 重映射正确。
  - gap 不满足 `min_gene_gap` 时抛出明确错误。
  - 无需拆分的染色体保留原区间，生成的 FASTA 序列全长一致。
- 覆盖边界：`split_size` 刚好等于染色体长度、`feature_type` 指定不存在的类型、FAI 中存在 0/负长度。 
