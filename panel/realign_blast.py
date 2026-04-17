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
from pandas.errors import EmptyDataError

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

LEGACY_BLAST_COLUMNS = [
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
BLAST_COLUMNS = [*LEGACY_BLAST_COLUMNS, "btop"]
BLAST_OUTFMT = f"6 {' '.join(BLAST_COLUMNS)}"
ALIGNMENT_COLUMNS = [
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
    "qstart",
    "qend",
    "sstart",
    "send",
    "btop",
]
MAPPING_COLUMNS = [
    *ALIGNMENT_COLUMNS,
    "match_ratio",
    "offset_fwd",
    "offset_rev",
    "alleles",
    "pos",
    "new_id",
    "pos_0",
    "n_count",
    "informative_len",
    "rank",
    "chr_map_status",
    "id_chrom_status",
    "selection_reason",
]
SELECTION_REPORT_COLUMNS = [
    "id",
    "new_id",
    "chrom",
    "pos",
    "strand",
    "alleles",
    "rank",
    "bitscore",
    "match_ratio",
    "align_len",
    "informative_len",
    "n_count",
    "mismatches",
    "gap_opens",
    "chr_map_status",
    "id_chrom_status",
    "selection_reason",
]
IUPAC_BASES = {
    "A": ("A",),
    "C": ("C",),
    "G": ("G",),
    "T": ("T",),
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("C", "G"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),
    "B": ("C", "G", "T"),
    "D": ("A", "G", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    "N": ("A", "T", "C", "G"),
}

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
        _run_cmd(f"bedtools slop -b {flank_size} -i {target_bed} -g {genome_fai} > {flank_bed}")
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
        _run_cmd(f"bedtools getfasta -fi {genome_fa} -fo {flank_fa} -bed {flank_bed} -nameOnly")
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
            f"blastn -query {query_fa} -db {db} "
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
    target_df = pd.read_table(target_bed, header=None, names=["target_start", "id"], usecols=[1, 3])
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
            hit_start(0-based), align_len, bitscore, mismatches, gap_opens,
            qstart, qend, sstart, send, btop
    """
    try:
        df = pd.read_table(blast_tsv, header=None)
    except EmptyDataError:
        return pd.DataFrame(columns=ALIGNMENT_COLUMNS)

    if df.empty:
        return pd.DataFrame(columns=ALIGNMENT_COLUMNS)

    if df.shape[1] == len(LEGACY_BLAST_COLUMNS):
        df.columns = LEGACY_BLAST_COLUMNS
        df["btop"] = pd.NA
    elif df.shape[1] == len(BLAST_COLUMNS):
        df.columns = BLAST_COLUMNS
    else:
        raise ValueError(
            f"BLAST 输出列数不符合预期: {blast_tsv}，实际 {df.shape[1]} 列，"
            f"预期 {len(LEGACY_BLAST_COLUMNS)} 或 {len(BLAST_COLUMNS)} 列"
        )

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

    return df.rename(
        columns={
            "qseqid": "id",
            "qlen": "query_len",
            "sseqid": "chrom",
            "length": "align_len",
            "mismatch": "mismatches",
            "gapopen": "gap_opens",
        }
    )[ALIGNMENT_COLUMNS].copy()


# ---------------------------------------------------------------------------
# Step 6: Infer target position
# ---------------------------------------------------------------------------


def _infer_position_without_gaps(row: pd.Series) -> int | None:
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


def _has_btop(value: object) -> bool:
    """判断 BLAST 结果是否包含可解析的 BTOP。"""
    return not pd.isna(value) and str(value).strip() != ""


def _tokenize_btop(btop: str) -> list[int | tuple[str, str]]:
    """解析 BLAST BTOP 字符串，返回 match 长度和 mismatch/gap pair。"""
    tokens: list[int | tuple[str, str]] = []
    index = 0
    while index < len(btop):
        if btop[index].isdigit():
            start = index
            while index < len(btop) and btop[index].isdigit():
                index += 1
            tokens.append(int(btop[start:index]))
            continue

        if index + 1 >= len(btop):
            raise ValueError(f"BTOP 格式错误，缺少成对碱基: {btop!r}")
        tokens.append((btop[index], btop[index + 1]))
        index += 2
    return tokens


def _offset_in_match_run(
    *,
    offset: int,
    q_current: int,
    q_step: int,
    run_length: int,
) -> bool:
    """判断 query offset 是否落在连续 match 区间内。"""
    if q_step > 0:
        return q_current <= offset < q_current + run_length
    return q_current >= offset > q_current - run_length


def _infer_position_from_btop(row: pd.Series) -> int | None:
    """
    根据 BTOP 精确推算目标位点坐标。

    BTOP 中 gap 位于目标位点之前时会校正 subject 坐标；目标位点本身
    落在 subject gap 上时返回 None，避免输出不可靠坐标。
    """
    btop = str(row["btop"]).strip()
    offset = int(row["offset_fwd"])
    qstart = int(row["qstart"]) - 1
    qend = int(row["qend"]) - 1
    sstart = int(row["sstart"]) - 1
    send = int(row["send"]) - 1

    q_step = 1 if qstart <= qend else -1
    s_step = 1 if sstart <= send else -1
    q_current = qstart
    s_current = sstart

    for token in _tokenize_btop(btop):
        if isinstance(token, int):
            if _offset_in_match_run(
                offset=offset,
                q_current=q_current,
                q_step=q_step,
                run_length=token,
            ):
                distance = (offset - q_current) * q_step
                return s_current + (distance * s_step) + 1
            q_current += token * q_step
            s_current += token * s_step
            continue

        query_base, subject_base = token
        consumes_query = query_base != "-"
        consumes_subject = subject_base != "-"

        if consumes_query and q_current == offset:
            if not consumes_subject:
                return None
            return s_current + 1

        if consumes_query:
            q_current += q_step
        if consumes_subject:
            s_current += s_step

    return None


def _infer_position(row: pd.Series) -> int | None:
    """推算目标位点在新基因组中的 1-based 坐标。"""
    if int(row["gap_opens"]) == 0 and not _has_btop(row.get("btop")):
        return _infer_position_without_gaps(row)
    if not _has_btop(row.get("btop")):
        raise ValueError("含 gap 的 BLAST 结果缺少 BTOP，无法精确推算坐标")
    return _infer_position_from_btop(row)


# ---------------------------------------------------------------------------
# Step 6b: Probe table parsing
# ---------------------------------------------------------------------------


def _probe_context(*, probe_table: Path | None, probe_id: object, flank: object) -> str:
    """生成探针表解析错误上下文。"""
    table = f"probe_table={probe_table}" if probe_table is not None else "probe_table=<unknown>"
    return f"{table}, id={probe_id!r}, Flank={flank!r}"


def _iupac_to_representative(sequence: str, *, context: str) -> str:
    """
    将 IUPAC 序列转换为代表性 ATGC 序列；保留 N 不变。

    BLAST 原生支持 query 中的 N（中性得分），为避免 “N→A” 造成的假阳性匹配，
      - N 字面写入序列；
      - R/Y/S/W/K/M/B/D/H/V 等其他 IUPAC 码继续按首位代表碱基展开，保持与原有行为一致。
    """
    converted: list[str] = []
    for base in sequence.upper():
        if base == "N":
            converted.append("N")
            continue
        try:
            converted.append(IUPAC_BASES[base][0])
        except KeyError as exc:
            raise ValueError(f"不支持的 IUPAC 碱基 {base!r}: {context}") from exc
    return "".join(converted)


def count_ns_in_fasta(fasta_path: Path) -> dict[str, int]:
    """
    统计 FASTA 中每条序列 N 的数量（大小写不敏感）。

    返回 {id: n_count}。id 取 FASTA header `>` 后的首个空格分隔符。
    """
    n_counts: dict[str, int] = {}
    current_id: str | None = None
    current_n = 0
    with open(fasta_path, encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    n_counts[current_id] = current_n
                current_id = line[1:].split()[0]
                current_n = 0
                continue
            current_n += line.upper().count("N")
    if current_id is not None:
        n_counts[current_id] = current_n
    return n_counts


def _validate_allele(allele: str, *, context: str) -> None:
    """验证等位基因仅包含 IUPAC 碱基或缺失标记。"""
    if allele == "":
        raise ValueError(f"等位基因不能为空: {context}")
    for base in allele.upper():
        if base != "-" and base not in IUPAC_BASES:
            raise ValueError(f"不支持的等位基因碱基 {base!r}: {context}")


def _representative_allele(alleles: list[str], *, context: str) -> str:
    """选择用于 BLAST query 的代表等位基因。"""
    for allele in alleles:
        if allele != "-":
            return _iupac_to_representative(allele, context=context)
    raise ValueError(f"等位基因不能全部为缺失标记: {context}")


def _parse_probe_flank(
    flank: object,
    *,
    probe_id: object = "<unknown>",
    probe_table: Path | None = None,
) -> tuple[str, list[tuple[int, str]]]:
    """
    从 Flank 字符串提取 query 序列、目标 offset 和 alleles 列表。

    支持 ``ACGT[A/G]TGCA``、``CGGCAA[Y]GACGCATTCG`` 和
    ``AC[A/G]TG[A/C]AAA`` 等含多个标记的 Flank。
    """
    context = _probe_context(probe_table=probe_table, probe_id=probe_id, flank=flank)
    if pd.isna(flank):
        raise ValueError(f"Flank 不能为空: {context}")

    flank_str = str(flank).strip()
    marker_matches = list(re.finditer(r"\[([^\[\]]+)\]", flank_str))
    if not marker_matches:
        raise ValueError(f"Flank 必须至少包含一个 [] 标记: {context}")

    sequence_parts: list[str] = []
    marker_rows: list[tuple[int, str]] = []
    cursor = 0
    sequence_length = 0

    for marker in marker_matches:
        literal_sequence = _iupac_to_representative(flank_str[cursor : marker.start()], context=context)
        sequence_parts.append(literal_sequence)
        sequence_length += len(literal_sequence)

        marker_text = marker.group(1).upper()
        if "/" in marker_text:
            alleles = marker_text.split("/")
            if len(alleles) != 2:
                raise ValueError(f"Flank 中 / 标记格式错误: {context}")
            for allele in alleles:
                _validate_allele(allele, context=context)
            center_sequence = _representative_allele(alleles, context=context)
            allele_text = "/".join(alleles)
        else:
            if len(marker_text) != 1 or marker_text not in IUPAC_BASES:
                raise ValueError(f"单碱基标记必须是一个有效 IUPAC 码: {context}")
            center_sequence = IUPAC_BASES[marker_text][0]
            allele_text = "/".join(IUPAC_BASES[marker_text])

        marker_rows.append((sequence_length, allele_text))
        sequence_parts.append(center_sequence)
        sequence_length += len(center_sequence)
        cursor = marker.end()

    trailing_sequence = _iupac_to_representative(flank_str[cursor:], context=context)
    sequence_parts.append(trailing_sequence)

    return "".join(sequence_parts), marker_rows


def _parse_probe_sequence(flank: str) -> str:
    """从 Flank 字符串提取纯序列（如 ``ACGT[A/G]TGCA`` → ``ACGTATGCA``）。"""
    return _parse_probe_flank(flank)[0]


def _parse_alleles(flank: str) -> str:
    """从 Flank 字符串提取 alleles（如 ``[A/G]`` → ``A/G``，``[Y]`` → ``C/T``）。"""
    return _parse_probe_flank(flank)[1][0][1]


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
    missing_columns = {"id", "Flank"} - set(df.columns)
    if missing_columns:
        raise ValueError(f"探针设计表缺少必需列 {sorted(missing_columns)}: {probe_table}")

    parsed_sequences = []
    offset_rows = []
    for row in df.itertuples():
        sequence, marker_rows = _parse_probe_flank(row.Flank, probe_id=row.id, probe_table=probe_table)
        parsed_sequences.append(sequence)
        seq_len = len(sequence)
        for offset_fwd, alleles in marker_rows:
            offset_rows.append(
                {
                    "id": row.id,
                    "offset_fwd": offset_fwd,
                    "offset_rev": seq_len - offset_fwd - 1,
                    "alleles": alleles,
                }
            )
    df["sequence"] = parsed_sequences

    probe_fasta = probe_table.with_suffix(".fa")
    with open(probe_fasta, "w", encoding="utf-8") as f:
        for row in df.itertuples():
            f.write(f">{row.id}\n{row.sequence}\n")

    return pd.DataFrame(offset_rows), probe_fasta


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


def load_id_chrom_map(id_chrom_map_file: Path) -> dict[str, set[str]]:
    """
    读取 ID 到目标染色体映射文件（2 列，无 header：id → 目标染色体）。

    返回 {id: {allowed_target_chrom, ...}} 字典。
    """
    try:
        df = pd.read_table(id_chrom_map_file, header=None, dtype=str)
    except EmptyDataError as exc:
        raise ValueError(f"ID 到染色体映射文件为空: {id_chrom_map_file}") from exc

    if df.empty:
        raise ValueError(f"ID 到染色体映射文件为空: {id_chrom_map_file}")
    if df.shape[1] != 2:
        raise ValueError(f"ID 到染色体映射文件必须为 2 列: {id_chrom_map_file}，实际 {df.shape[1]} 列")

    df.columns = ["id", "chrom"]
    df["id"] = df["id"].str.strip()
    df["chrom"] = df["chrom"].str.strip()
    invalid_rows = df["id"].isna() | df["chrom"].isna() | (df["id"] == "") | (df["chrom"] == "")
    if invalid_rows.any():
        row_numbers = ", ".join(str(index + 1) for index in df.index[invalid_rows].tolist())
        raise ValueError(f"ID 到染色体映射文件存在空 id/chrom: {id_chrom_map_file}，行: {row_numbers}")

    return df.groupby("id")["chrom"].apply(set).to_dict()


def load_source_chroms(target_bed: Path) -> pd.DataFrame:
    """从 target BED 读取每个 id 的源染色体。"""
    return pd.read_table(target_bed, header=None, names=["source_chrom", "id"], usecols=[0, 3])


def _empty_mapping_df() -> pd.DataFrame:
    """返回带输出所需列的空映射表。"""
    return pd.DataFrame(columns=MAPPING_COLUMNS)


def _filter_by_id_chrom_map(
    df: pd.DataFrame,
    id_chrom_map: dict[str, set[str]],
) -> pd.DataFrame:
    """
    按 id -> 目标染色体映射过滤 BLAST hit。

    仅对 id_chrom_map 中登记的 id 施加严格过滤（只保留 chrom 位于允许集合的 hit）；
    未登记的 id 整体放行，后续走默认 best-hit 排序策略。对缺失登记的 id 会输出
    warning 日志，避免静默丢数据。
    """
    if df.empty:
        return df.copy()

    ids_in_df = set(df["id"].astype(str).unique())
    unmapped_ids = ids_in_df - id_chrom_map.keys()
    if unmapped_ids:
        sample = sorted(unmapped_ids)[:5]
        suffix = ", ..." if len(unmapped_ids) > len(sample) else ""
        logger.warning(
            f"{len(unmapped_ids)} / {len(ids_in_df)} 个 id 未在 id_chrom_map 中登记，"
            f"将按默认 best-hit 策略保留其所有 hit（示例: {sample}{suffix}）"
        )

    def _is_allowed(row: pd.Series) -> bool:
        allowed = id_chrom_map.get(str(row["id"]))
        if allowed is None:
            return True
        return str(row["chrom"]) in allowed

    has_allowed_chrom = df.apply(_is_allowed, axis=1)
    return df[has_allowed_chrom].copy()


def _ensure_btop_for_gaps(df: pd.DataFrame, *, max_gap_opens: int) -> None:
    """允许 gap 时，确保含 gap 的 hit 都有 BTOP 可用于精确坐标推断。"""
    if max_gap_opens <= 0 or df.empty:
        return
    gap_rows = df["gap_opens"] > 0
    if gap_rows.any() and "btop" not in df.columns:
        raise ValueError("允许 gap 时需要 BLAST BTOP 列；请使用 --force 重新生成 BLAST 结果")
    missing_btop = gap_rows & ~df["btop"].map(_has_btop)
    if missing_btop.any():
        raise ValueError("允许 gap 时需要 BLAST BTOP 列；请使用 --force 重新生成 BLAST 结果")


def build_id_mapping(
    alignments: pd.DataFrame,
    offsets: pd.DataFrame,
    *,
    match_ratio_cutoff: float = 0.9,
    max_hits: int = 1,
    chr_map: dict[str, set[str]] | None = None,
    source_chroms: pd.DataFrame | None = None,
    id_chrom_map: dict[str, set[str]] | None = None,
    max_gap_opens: int = 0,
    n_counts: dict[str, int] | None = None,
) -> pd.DataFrame:
    """
    筛选最佳比对、计算新坐标、生成 ID 映射表。

    筛选策略：
      1. 过滤 gap_open 数量超过 max_gap_opens 的比对
      2. match_ratio = align_len / informative_len，其中 informative_len = query_len - n_count；
         滤掉 match_ratio <= match_ratio_cutoff 的比对。informative_len ≤ 0（全 N）的 id 被丢弃并 warn。
      3. 若提供 id_chrom_map，仅对登记过的 id 严格过滤到允许染色体上的 hit；
         未登记的 id 整体放行走默认排序策略，并输出 warning 日志
      4. 若提供 chr_map，优先保留目标染色体匹配的 hit；无匹配时回退到最佳 hit
      5. 每个 id 按 bitscore 降序 → mismatches 升序排列
      6. 每个 id 最多保留 max_hits 条

    n_counts: {id: N 碱基数量}。未提供时视为 0，等价于旧行为。
    """
    if max_gap_opens < 0:
        raise ValueError(f"max_gap_opens 不能小于 0: {max_gap_opens}")

    df = alignments.copy()

    # 先过滤：仅保留 gap 数量达标的 hit
    df = df[df["gap_opens"] <= max_gap_opens].copy()
    _ensure_btop_for_gaps(df, max_gap_opens=max_gap_opens)

    # 计算 informative_len（排除 N 后的有效 query 长度）作为 match_ratio 分母
    n_counts = n_counts or {}
    df["n_count"] = df["id"].astype(str).map(n_counts).fillna(0).astype(int)
    df["informative_len"] = df["query_len"] - df["n_count"]

    too_short = df["informative_len"] <= 0
    if too_short.any():
        bad_ids = sorted(df.loc[too_short, "id"].astype(str).unique())
        sample = bad_ids[:5]
        suffix = ", ..." if len(bad_ids) > len(sample) else ""
        logger.warning(
            f"{len(bad_ids)} 个 id 的 query 信息量为 0（全为 N），已丢弃其所有 hit，"
            f"示例: {sample}{suffix}"
        )
        df = df[~too_short].copy()

    if df.empty:
        return _empty_mapping_df()

    df["match_ratio"] = df["align_len"] / df["informative_len"]
    df = df[df["match_ratio"] > match_ratio_cutoff].copy()

    if df.empty:
        return _empty_mapping_df()

    if id_chrom_map is not None:
        df = _filter_by_id_chrom_map(df, id_chrom_map)
        if df.empty:
            return _empty_mapping_df()

    # 注入 id_chrom_map 状态（供报告使用）
    if id_chrom_map is not None:
        df["id_chrom_status"] = df["id"].astype(str).apply(
            lambda i: "strict" if i in id_chrom_map else "fallback"
        )
    else:
        df["id_chrom_status"] = "n/a"

    # chr_map 优先排序：匹配的 hit 排在前面
    if chr_map is not None and source_chroms is not None:
        df = df.merge(source_chroms, on="id", how="left")
        df["chr_matched"] = df.apply(
            lambda r: r["chrom"] in chr_map.get(str(r["source_chrom"]), set()),
            axis=1,
        )
        df["chr_map_status"] = df["chr_matched"].map({True: "matched", False: "fallback"})
        # chr_matched=True 排前面（ascending=True 时 False<True，所以用 ascending=False）
        df = df.sort_values(
            ["id", "chr_matched", "bitscore", "mismatches"],
            ascending=[True, False, False, True],
            kind="mergesort",
        )
        df = df.drop(columns=["source_chrom", "chr_matched"])
    else:
        df["chr_map_status"] = "n/a"
        df = df.sort_values(
            ["id", "bitscore", "mismatches"],
            ascending=[True, False, True],
            kind="mergesort",
        )

    # 排序后分配 rank（1 = 当前 id 下的最优 hit）
    df["rank"] = df.groupby("id", sort=False).cumcount() + 1
    df = df.groupby("id", sort=False).head(max_hits)

    # 合并偏移量并推算新坐标
    df = df.merge(offsets, on="id")
    df["pos"] = df.apply(_infer_position, axis=1)
    df = df.dropna(subset=["pos"])
    if df.empty:
        return _empty_mapping_df()
    df["pos"] = df["pos"].astype(int)
    df["new_id"] = df["chrom"].astype(str) + "_" + df["pos"].astype(str)
    df["pos_0"] = df["pos"] - 1
    df["selection_reason"] = df.apply(_compose_selection_reason, axis=1)

    return df


def _compose_selection_reason(row: pd.Series) -> str:
    """生成单条待选中 hit 的人类可读选择理由。"""
    parts = [
        f"rank #{int(row['rank'])}",
        f"bitscore={float(row['bitscore']):.1f}",
        f"match_ratio={float(row['match_ratio']):.3f}",
        f"mismatches={int(row['mismatches'])}",
        f"gap_opens={int(row['gap_opens'])}",
    ]
    chr_status = row.get("chr_map_status", "n/a")
    if chr_status == "matched":
        parts.append("chr_map: 源→目标染色体匹配")
    elif chr_status == "fallback":
        parts.append("chr_map: 未匹配，取跨染色体最优")
    id_status = row.get("id_chrom_status", "n/a")
    if id_status == "strict":
        parts.append("id_chrom_map: 严格限制通过")
    elif id_status == "fallback":
        parts.append("id_chrom_map: 未登记，走默认 best-hit")
    return "; ".join(parts)


# ---------------------------------------------------------------------------
# Step 8: Write outputs
# ---------------------------------------------------------------------------


def write_outputs(mapping_df: pd.DataFrame, out_prefix: Path) -> None:
    """写出 idmap.tsv / target.bed / pos.tsv / selection.tsv 四个结果文件。"""
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

    # selection.tsv：对外可解释的选择报告
    write_selection_report(mapping_df, idmap_path)


def write_selection_report(mapping_df: pd.DataFrame, out_prefix: Path) -> Path:
    """
    写出 *.selection.tsv 报告：每条选中 hit 的完整数值指标与选择理由。

    输出列见 ``SELECTION_REPORT_COLUMNS``，含表头；缺失 alleles 时默认填 ``-/-``。
    返回写出的报告路径。
    """
    report_path = out_prefix.with_suffix(".selection.tsv")
    df = mapping_df.copy()
    if "alleles" not in df.columns:
        df["alleles"] = "-/-"
    df.to_csv(
        report_path,
        sep="\t",
        index=False,
        columns=SELECTION_REPORT_COLUMNS,
    )
    logger.info(f"选择依据报告: {report_path}（{len(df)} 条记录）")
    return report_path


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

    # 5. 计算偏移量和 N 统计
    offsets = compute_offsets(target_bed, flank_bed)
    n_counts = count_ns_in_fasta(flank_fa)

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
        n_counts=n_counts,
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
          --max-hits 1 --match-ratio-cutoff 0.95 \\
          --id-chrom-map id_chrom.tsv --max-gap-opens 2
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
    id_chrom_map_file: Annotated[
        Path | None,
        typer.Option(
            "--id-chrom-map",
            help="ID 到目标染色体映射文件（2 列：id → 目标染色体）；仅对登记过的 id 严格过滤，未登记的 id 走默认 best-hit 策略",
        ),
    ] = None,
    max_gap_opens: Annotated[int, typer.Option(help="允许的最大 gap opening 数量；大于 0 时需要 BTOP 列")] = 2,
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
    n_counts = count_ns_in_fasta(probe_fa)
    id_chrom_map = None
    if id_chrom_map_file is not None:
        logger.info(f"加载 ID 到染色体映射: {id_chrom_map_file}")
        id_chrom_map = load_id_chrom_map(id_chrom_map_file)

    mapping_df = build_id_mapping(
        alignments,
        offsets,
        match_ratio_cutoff=match_ratio_cutoff,
        max_hits=max_hits,
        id_chrom_map=id_chrom_map,
        max_gap_opens=max_gap_opens,
        n_counts=n_counts,
    )

    # 4. 输出结果
    write_outputs(mapping_df, blast_tsv)


if __name__ == "__main__":
    app()
