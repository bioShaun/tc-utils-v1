"""
基于 BLAST 比对将目标位点坐标从源参考基因组映射到目标参考基因组。

\b
子命令:
  from-bed          从 BED/VCF 坐标文件 + 源基因组提取侧翼序列，重定位到目标基因组
  from-probe-table  从探针设计表（含 Flank 序列）重定位到目标基因组

\b
输出格式:
  - *.idmap.tsv:      old_id \\t new_id（无 header）
  - *.target.bed:     BED6（chrom, start, end, id, bitscore, strand）
  - *.pos.tsv:        chrom \\t pos \\t alleles \\t id（有 header）
"""

from __future__ import annotations

import re
import subprocess
from enum import Enum
from inspect import cleandoc
from pathlib import Path
from typing import Annotated

import pandas as pd
import typer
from loguru import logger

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BLAST_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
]
BLAST_OUTFMT = f"6 {' '.join(BLAST_COLUMNS)}"

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

app = typer.Typer(
    help="基于 BLAST 比对将目标位点坐标从源参考基因组映射到目标参考基因组。",
    no_args_is_help=True,
)


class TargetType(str, Enum):
    bed = "bed"
    vcf = "vcf"


# ---------------------------------------------------------------------------
# Step 1: BED from VCF
# ---------------------------------------------------------------------------


def bed_from_vcf(vcf: Path) -> Path:
    """将 VCF 转换为 4 列 BED（chrom, start, end, id）。"""
    df = pd.read_table(
        vcf,
        header=None,
        usecols=[0, 1],
        comment="#",
        names=["chrom", "end"],
    )
    df["start"] = df["end"] - 1
    df["id"] = df["chrom"].astype(str) + "_" + df["end"].astype(str)
    bed_path = vcf.with_suffix(".bed")
    df.to_csv(
        bed_path,
        sep="\t",
        index=False,
        header=False,
        columns=["chrom", "start", "end", "id"],
    )
    return bed_path


# ---------------------------------------------------------------------------
# Step 2: Flank BED / FASTA
# ---------------------------------------------------------------------------


def _run_cmd(cmd: str) -> None:
    """执行 shell 命令，失败时抛出异常并包含 stderr 信息。"""
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"命令执行失败: {cmd}\nstderr: {result.stderr.strip()}")


def make_flank_bed(
    target_bed: Path,
    flank_size: int,
    genome_fai: Path,
    *,
    force: bool = False,
) -> Path:
    """用 bedtools slop 扩展目标位点侧翼区域。"""
    flank_bed = target_bed.with_suffix(f".flank{flank_size}.bed")
    if force or not flank_bed.is_file():
        _run_cmd(
            f"bedtools slop -b {flank_size} -i {target_bed} -g {genome_fai} > {flank_bed}"
        )
    return flank_bed


def extract_flank_fasta(
    flank_bed: Path,
    genome_fa: Path,
    *,
    force: bool = False,
) -> Path:
    """用 bedtools getfasta 提取侧翼序列。"""
    flank_fa = flank_bed.with_suffix(".fa")
    if force or not flank_fa.is_file():
        _run_cmd(
            f"bedtools getfasta -fi {genome_fa} -fo {flank_fa} -bed {flank_bed} -nameOnly"
        )
    return flank_fa


# ---------------------------------------------------------------------------
# Step 3: BLAST
# ---------------------------------------------------------------------------


def run_blastn(
    query_fa: Path,
    db: Path,
    *,
    threads: int = 16,
    evalue: float = 1e-10,
    max_target_seqs: int = 10,
    force: bool = False,
) -> Path:
    """执行 blastn 比对，输出 tabular 格式。"""
    blast_tsv = query_fa.with_suffix(".blast.tsv")
    if force or not blast_tsv.is_file():
        _run_cmd(
            f'blastn -query {query_fa} -db {db} '
            f'-outfmt "{BLAST_OUTFMT}" '
            f"-num_threads {threads} "
            f"-evalue {evalue} "
            f"-max_target_seqs {max_target_seqs} "
            f"-out {blast_tsv}"
        )
    return blast_tsv


# ---------------------------------------------------------------------------
# Step 4: Offset computation
# ---------------------------------------------------------------------------


def compute_offsets(target_bed: Path, flank_bed: Path) -> pd.DataFrame:
    """
    计算目标位点在侧翼序列中的前向/反向偏移量。

    offset_fwd: 目标位点距 flank 起始的碱基数（正链时使用）
    offset_rev: 目标位点距 flank 末尾的碱基数（负链时使用）
    """
    target_df = pd.read_table(
        target_bed, header=None, names=["target_start", "id"], usecols=[1, 3]
    )
    flank_df = pd.read_table(
        flank_bed,
        header=None,
        names=["flank_start", "flank_end", "id"],
        usecols=[1, 2, 3],
    )
    merged = target_df.merge(flank_df, on="id")
    merged["offset_fwd"] = merged["target_start"] - merged["flank_start"]
    merged["offset_rev"] = merged["flank_end"] - merged["target_start"] - 1
    return merged[["id", "offset_fwd", "offset_rev"]].copy()


# ---------------------------------------------------------------------------
# Step 5: Parse BLAST → alignment DataFrame
# ---------------------------------------------------------------------------


def parse_blast_results(blast_tsv: Path) -> pd.DataFrame:
    """
    读取 BLAST tabular 输出，标准化正负链坐标。

    返回列: id, query_len, query_start(0-based), strand, chrom,
            hit_start(0-based), align_len, bitscore, mismatches, gap_opens
    """
    df = pd.read_table(blast_tsv, header=None, names=BLAST_COLUMNS)

    # 判断链方向：sstart < send 为正链
    is_positive = df["sstart"] < df["send"]

    # 正链: query_start = qstart - 1 (BLAST 1-based → 0-based)
    #        hit_start   = sstart - 1
    # 负链: query_start = qlen - qend (反向互补时从 query 尾部算起)
    #        hit_start   = send - 1   (交换后较小值)
    df["strand"] = "+"
    df.loc[~is_positive, "strand"] = "-"

    df["query_start"] = df["qstart"] - 1
    df["hit_start"] = df["sstart"] - 1
    df.loc[~is_positive, "query_start"] = df.loc[~is_positive, "qlen"] - df.loc[~is_positive, "qend"]
    df.loc[~is_positive, "hit_start"] = df.loc[~is_positive, "send"] - 1

    return df.rename(columns={
        "qseqid": "id",
        "qlen": "query_len",
        "sseqid": "chrom",
        "length": "align_len",
        "mismatch": "mismatches",
        "gapopen": "gap_opens",
    })[
        [
            "id",
            "query_len",
            "query_start",
            "strand",
            "chrom",
            "hit_start",
            "align_len",
            "bitscore",
            "mismatches",
            "gap_opens",
        ]
    ].copy()


# ---------------------------------------------------------------------------
# Step 6: Infer target position
# ---------------------------------------------------------------------------


def _infer_position(row: pd.Series) -> int | None:
    """
    根据比对位置和偏移量推算目标位点在新基因组的 1-based 坐标。

    仅处理无 gap 的比对（gap_opens == 0），有 gap 时返回 None。

    坐标推算逻辑（以正链为例）:
      - hit_start (0-based) 是侧翼序列在目标基因组的起始比对位置
      - query_start (0-based) 是侧翼序列中实际参与比对的起始偏移
      - offset_fwd 是目标位点在完整侧翼序列中的偏移
      - 新坐标 = hit_start + (offset_fwd - query_start) + 1
    """
    if row["gap_opens"] > 0:
        return None
    offset = row["offset_fwd"] if row["strand"] == "+" else row["offset_rev"]
    # 检查 offset 是否在比对覆盖范围内
    if offset < row["query_start"]:
        return None
    if offset >= row["query_start"] + row["align_len"]:
        return None
    return row["hit_start"] + (offset - row["query_start"]) + 1


# ---------------------------------------------------------------------------
# Step 6b: Probe table parsing
# ---------------------------------------------------------------------------


def _parse_probe_sequence(flank: str) -> str:
    """从 Flank 字符串提取纯序列（如 ``ACGT[A/G]TGCA`` → ``ACGTATGCA``）。"""
    left = flank.split("[")[0]
    center = flank.split("[")[1][0]
    right = flank.split("]")[1]
    return f"{left}{center}{right}"


def _parse_alleles(flank: str) -> str:
    """从 Flank 字符串提取 alleles（如 ``[A/G]`` → ``A/G``）。"""
    match = re.search(r"\[([ACGT\-]+)/([ACGT\-]+)\]", flank)
    if match:
        return f"{match.group(1)}/{match.group(2)}"
    return "-/-"


def fasta_and_offsets_from_probe_table(
    probe_table: Path,
) -> tuple[pd.DataFrame, Path]:
    """
    从探针设计表生成 FASTA 和偏移量表。

    探针设计表需包含 ``id`` 和 ``Flank`` 列，Flank 格式如 ``ACGT[A/G]TGCA``。

    返回:
      - offsets DataFrame（列：id, offset_fwd, offset_rev, alleles）
      - probe FASTA 文件路径
    """
    df = pd.read_table(probe_table)
    df["sequence"] = df["Flank"].map(_parse_probe_sequence)

    probe_fasta = probe_table.with_suffix(".fa")
    with open(probe_fasta, "w", encoding="utf-8") as f:
        for row in df.itertuples():
            f.write(f">{row.id}\n{row.sequence}\n")

    seq_len = df["sequence"].str.len()
    df["offset_fwd"] = df["Flank"].map(lambda x: x.index("["))
    df["offset_rev"] = seq_len - df["offset_fwd"] - 1
    df["alleles"] = df["Flank"].map(_parse_alleles)

    return df[["id", "offset_fwd", "offset_rev", "alleles"]].copy(), probe_fasta


# ---------------------------------------------------------------------------
# Step 7: Build ID mapping
# ---------------------------------------------------------------------------


def load_chr_map(chr_map_file: Path) -> dict[str, set[str]]:
    """
    读取染色体映射文件（2 列，无 header：源染色体 → 目标染色体）。

    返回 {source_chrom: {allowed_target_chrom, ...}} 字典。
    同一源染色体可映射到多个目标染色体。
    """
    df = pd.read_table(chr_map_file, header=None, names=["source", "target"])
    return df.groupby("source")["target"].apply(set).to_dict()


def load_source_chroms(target_bed: Path) -> pd.DataFrame:
    """从 target BED 读取每个 id 的源染色体。"""
    return pd.read_table(
        target_bed, header=None, names=["source_chrom", "id"], usecols=[0, 3]
    )


def build_id_mapping(
    alignments: pd.DataFrame,
    offsets: pd.DataFrame,
    *,
    match_ratio_cutoff: float = 0.9,
    max_hits: int = 1,
    chr_map: dict[str, set[str]] | None = None,
    source_chroms: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """
    筛选最佳比对、计算新坐标、生成 ID 映射表。

    筛选策略：
      1. 过滤含 gap 的比对（gap_opens > 0）
      2. 过滤 align_len / query_len <= match_ratio_cutoff 的比对
      3. 若提供 chr_map，优先保留目标染色体匹配的 hit；无匹配时回退到最佳 hit
      4. 每个 id 按 bitscore 降序 → mismatches 升序排列
      5. 每个 id 最多保留 max_hits 条
    """
    df = alignments.copy()

    # 先过滤：仅保留无 gap 且比对比率达标的 hit
    df = df[df["gap_opens"] == 0].copy()
    df["match_ratio"] = df["align_len"] / df["query_len"]
    df = df[df["match_ratio"] > match_ratio_cutoff].copy()

    if df.empty:
        return df

    # chr_map 优先排序：匹配的 hit 排在前面
    if chr_map is not None and source_chroms is not None:
        df = df.merge(source_chroms, on="id", how="left")
        df["chr_matched"] = df.apply(
            lambda r: r["chrom"] in chr_map.get(str(r["source_chrom"]), set()),
            axis=1,
        )
        # chr_matched=True 排前面（ascending=True 时 False<True，所以用 ascending=False）
        df = df.sort_values(
            ["id", "chr_matched", "bitscore", "mismatches"],
            ascending=[True, False, False, True],
        )
        df = df.drop(columns=["source_chrom", "chr_matched"])
    else:
        df = df.sort_values(
            ["id", "bitscore", "mismatches"], ascending=[True, False, True]
        )

    df = df.groupby("id", sort=False).head(max_hits)

    # 合并偏移量并推算新坐标
    df = df.merge(offsets, on="id")
    df["pos"] = df.apply(_infer_position, axis=1)
    df = df.dropna(subset=["pos"])
    df["pos"] = df["pos"].astype(int)
    df["new_id"] = df["chrom"].astype(str) + "_" + df["pos"].astype(str)
    df["pos_0"] = df["pos"] - 1

    return df


# ---------------------------------------------------------------------------
# Step 8: Write outputs
# ---------------------------------------------------------------------------


def write_outputs(mapping_df: pd.DataFrame, out_prefix: Path) -> None:
    """写出 idmap.tsv / target.bed / pos.tsv 三个结果文件。"""
    # idmap.tsv
    idmap_path = out_prefix.with_suffix(".idmap.tsv")
    mapping_df.to_csv(
        idmap_path,
        sep="\t",
        header=False,
        index=False,
        columns=["id", "new_id"],
    )
    logger.info(f"ID 映射表: {idmap_path}（{len(mapping_df)} 条记录）")

    # target.bed (sorted)
    bed_path = idmap_path.with_suffix(".target.bed")
    sorted_df = mapping_df.sort_values(["chrom", "pos"])
    sorted_df.to_csv(
        bed_path,
        sep="\t",
        index=False,
        header=False,
        columns=["chrom", "pos_0", "pos", "id", "bitscore", "strand"],
    )
    logger.info(f"目标 BED: {bed_path}")

    # pos.tsv（含 alleles 列，若有）
    pos_path = idmap_path.with_suffix(".pos.tsv")
    if "alleles" not in mapping_df.columns:
        mapping_df = mapping_df.copy()
        mapping_df["alleles"] = "-/-"
    mapping_df.to_csv(
        pos_path,
        sep="\t",
        index=False,
        columns=["chrom", "pos", "alleles", "id"],
    )
    logger.info(f"位点坐标: {pos_path}")


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

FROM_BED_HELP = cleandoc(
    """
    从 BED/VCF 坐标文件 + 源基因组提取侧翼序列，BLAST 比对到目标基因组，推算新坐标。

    \b
    使用示例:
      python panel/realign_blast.py from-bed target.bed genome.fa blast_db

      python panel/realign_blast.py from-bed target.vcf genome.fa blast_db \\
          --target-type vcf --flank-size 150
    """
)


@app.command(help=FROM_BED_HELP)
def from_bed(
    target_file: Annotated[Path, typer.Argument(help="目标位点文件（BED 或 VCF）")],
    genome_fa: Annotated[Path, typer.Argument(help="源参考基因组 FASTA（用于提取侧翼序列）")],
    blast_db: Annotated[Path, typer.Argument(help="目标基因组 BLAST 数据库前缀")],
    flank_size: Annotated[int, typer.Option(help="侧翼区域大小（bp）")] = 100,
    threads: Annotated[int, typer.Option(help="BLAST 并行线程数")] = 16,
    evalue: Annotated[float, typer.Option(help="BLAST E-value 阈值")] = 1e-10,
    match_ratio_cutoff: Annotated[float, typer.Option(help="比对长度/序列长度最低比率")] = 0.9,
    max_target_seqs: Annotated[int, typer.Option(help="BLAST 每条 query 最大目标序列数")] = 10,
    target_type: Annotated[TargetType, typer.Option(help="输入文件类型")] = TargetType.bed,
    force: Annotated[bool, typer.Option(help="强制重新运行所有中间步骤")] = False,
    max_hits: Annotated[int, typer.Option(help="每个位点最多保留的 best hit 数量")] = 3,
    chr_map_file: Annotated[
        Path | None,
        typer.Option("--chr-map", help="染色体映射文件（2 列：源染色体 → 目标染色体），优先保留匹配的 hit"),
    ] = None,
) -> None:
    genome_fai = genome_fa.parent / f"{genome_fa.name}.fai"

    # 1. 确定 target BED
    if target_type == TargetType.vcf:
        logger.info("从 VCF 生成 BED ...")
        target_bed = bed_from_vcf(target_file)
    else:
        target_bed = target_file

    # 2. 扩展侧翼区域
    logger.info(f"生成 ±{flank_size}bp 侧翼 BED ...")
    flank_bed = make_flank_bed(target_bed, flank_size, genome_fai, force=force)

    # 3. 提取侧翼 FASTA
    logger.info("提取侧翼序列 FASTA ...")
    flank_fa = extract_flank_fasta(flank_bed, genome_fa, force=force)

    # 4. BLAST 比对
    logger.info("运行 blastn 比对 ...")
    blast_tsv = run_blastn(
        flank_fa,
        blast_db,
        threads=threads,
        evalue=evalue,
        max_target_seqs=max_target_seqs,
        force=force,
    )

    # 5. 计算偏移量
    offsets = compute_offsets(target_bed, flank_bed)

    # 6. 加载 chr-map（如果提供）
    chr_map = None
    source_chroms = None
    if chr_map_file is not None:
        logger.info(f"加载染色体映射: {chr_map_file}")
        chr_map = load_chr_map(chr_map_file)
        source_chroms = load_source_chroms(target_bed)

    # 7. 解析 BLAST 结果并生成映射
    logger.info("解析比对结果、推算新坐标 ...")
    alignments = parse_blast_results(blast_tsv)
    mapping_df = build_id_mapping(
        alignments,
        offsets,
        match_ratio_cutoff=match_ratio_cutoff,
        max_hits=max_hits,
        chr_map=chr_map,
        source_chroms=source_chroms,
    )

    # 8. 输出结果
    write_outputs(mapping_df, blast_tsv)


FROM_PROBE_TABLE_HELP = cleandoc(
    """
    从探针设计表（含 Flank 序列）提取探针序列，BLAST 比对到目标基因组，推算新坐标。

    探针设计表需包含 id 和 Flank 列，Flank 格式如 ACGT[A/G]TGCA。

    \b
    使用示例:
      python panel/realign_blast.py from-probe-table probe.tsv blast_db

      python panel/realign_blast.py from-probe-table probe.tsv blast_db \\
          --max-hits 1 --match-ratio-cutoff 0.95
    """
)


@app.command(help=FROM_PROBE_TABLE_HELP)
def from_probe_table(
    probe_table: Annotated[Path, typer.Argument(help="探针设计表（含 id 和 Flank 列）")],
    blast_db: Annotated[Path, typer.Argument(help="目标基因组 BLAST 数据库前缀")],
    threads: Annotated[int, typer.Option(help="BLAST 并行线程数")] = 16,
    evalue: Annotated[float, typer.Option(help="BLAST E-value 阈值")] = 1e-10,
    match_ratio_cutoff: Annotated[float, typer.Option(help="比对长度/序列长度最低比率")] = 0.9,
    max_target_seqs: Annotated[int, typer.Option(help="BLAST 每条 query 最大目标序列数")] = 10,
    force: Annotated[bool, typer.Option(help="强制重新运行 BLAST")] = False,
    max_hits: Annotated[int, typer.Option(help="每个位点最多保留的 best hit 数量")] = 3,
) -> None:
    # 1. 从探针设计表生成 FASTA 和偏移量
    logger.info("解析探针设计表，生成 FASTA ...")
    offsets, probe_fa = fasta_and_offsets_from_probe_table(probe_table)

    # 2. BLAST 比对
    logger.info("运行 blastn 比对 ...")
    blast_tsv = run_blastn(
        probe_fa,
        blast_db,
        threads=threads,
        evalue=evalue,
        max_target_seqs=max_target_seqs,
        force=force,
    )

    # 3. 解析 BLAST 结果并生成映射
    logger.info("解析比对结果、推算新坐标 ...")
    alignments = parse_blast_results(blast_tsv)
    mapping_df = build_id_mapping(
        alignments,
        offsets,
        match_ratio_cutoff=match_ratio_cutoff,
        max_hits=max_hits,
    )

    # 4. 输出结果
    write_outputs(mapping_df, blast_tsv)


if __name__ == "__main__":
    app()
