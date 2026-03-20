"""
根据设计表生成 panel BED 文件，并可按 split.bed 拆分坐标。

使用示例:
    # 默认输出
    python chip/panel2bed.py design.tsv genome.fa.fai panel_v1 out_dir

    # 输出 split 版本 BED（按 split.genome.fa.fai 排序）
    python chip/panel2bed.py design.tsv genome.fa.fai panel_v1 out_dir --split-bed split.bed --split-genome-fai split.genome.fa.fai
"""

from pathlib import Path
from inspect import cleandoc
from typing import Annotated, Iterable, Sequence, Tuple

import pandas as pd
import typer

REQUIRED_COLUMNS = ("chrom", "pos")
MODULE_HELP = cleandoc(
    """
    根据设计表生成 panel BED 文件，并可按 split.bed 拆分坐标。

    \b
    使用示例:
    \b
      1) 默认输出
      python chip/panel2bed.py \\
        design.tsv \\
        genome.fa.fai \\
        panel_v1 \\
        out_dir

    \b
      2) 输出 split 版本 BED（按 split.genome.fa.fai 排序）
      python chip/panel2bed.py \\
        design.tsv \\
        genome.fa.fai \\
        panel_v1 \\
        out_dir \\
        --split-bed split.bed \\
        --split-genome-fai split.genome.fa.fai

    \b
    输出格式:
      - <probe_id>.id: 1 列，pos_id（chrom_pos）
      - <probe_id>.bed: 3 列（chrom, start, end），未启用 split 时输出
      - <probe_id>.snpcalling.bed: 3 列（chrom, start, end），未启用 split 时输出
      - <probe_id>.bed: 3 列（new_chrom, new_start, new_end），启用 split 时输出
      - <probe_id>.snpcalling.bed: 3 列（new_chrom, new_start, new_end），启用 split 时输出
    """
)
app = typer.Typer(help=MODULE_HELP, no_args_is_help=True)


def merge_intervals(group: pd.DataFrame, start_col: str, end_col: str) -> pd.DataFrame:
    """
    Merge overlapping [start_col, end_col] intervals. Assumes closed intervals.
    """
    if group.empty:
        return group[[start_col, end_col]]

    ordered = group.sort_values(start_col).reset_index(drop=True)
    ordered[[start_col, end_col]] = ordered[[start_col, end_col]].apply(
        pd.to_numeric, errors="coerce"
    )
    merged = []
    current_start = ordered.loc[0, start_col]
    current_end = ordered.loc[0, end_col]

    for _, row in ordered.iloc[1:].iterrows():
        start = row[start_col]
        end = row[end_col]
        if pd.isna(start) or pd.isna(end):
            continue
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end
    merged.append((current_start, current_end))

    return pd.DataFrame(merged, columns=[start_col, end_col])


def get_flank_start_end(
    probe_start: int, probe_end: int, flank_size: int, chr_size: int
) -> Tuple[int, int]:
    """
    Expand probe interval to a requested flank size while respecting chromosome bounds.
    """
    probe_len = probe_end - probe_start
    extend_len = (flank_size - probe_len) // 2
    flank_start = max(0, probe_start - extend_len)
    flank_end = min(chr_size, probe_end + extend_len)
    return flank_start, flank_end


def ensure_columns(df: pd.DataFrame, required: Iterable[str]) -> None:
    """
    Ensure that the design table has the required columns.
    """
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"设计表缺少字段: {', '.join(missing)}")


def load_chrom_sizes(genome_fai: Path) -> pd.DataFrame:
    """Load chromosome sizes from FASTA index."""
    chrom_df = pd.read_table(
        genome_fai, header=None, names=["chrom", "chrom_size"], usecols=[0, 1]
    )
    chrom_df["chrom"] = chrom_df["chrom"].astype(str)
    chrom_df["chrom_size"] = pd.to_numeric(
        chrom_df["chrom_size"], errors="raise"
    ).astype(int)
    return chrom_df


def load_design_table(design_table: Path) -> pd.DataFrame:
    """Load design table."""
    if not design_table.exists():
        raise FileNotFoundError(f"找不到设计表: {design_table}")

    suffix = design_table.suffix.lower()
    if suffix == ".tsv":
        df = pd.read_csv(design_table, sep="\t")
    elif suffix in {".xlsx", ".xls"}:
        df = pd.read_excel(design_table)
    else:
        raise ValueError("设计表只支持 .tsv / .xlsx / .xls 格式")
    return df


def prepare_probe_dataframe(df: pd.DataFrame, chrom_df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize chromosome labels, add pos_id、pos_0 并按照染色体顺序排序。
    """
    ensure_columns(df, REQUIRED_COLUMNS)
    df = df.copy()
    df["chrom"] = df["chrom"].astype(str)
    df["pos"] = pd.to_numeric(df["pos"], errors="raise").astype(int)

    # Bring in chromosome size to bound default probe intervals.
    df = df.merge(chrom_df, how="left", on="chrom")
    if df["chrom_size"].isna().any():
        invalid = df.loc[df["chrom_size"].isna(), "chrom"].unique()
        raise ValueError(f"缺少染色体大小: {', '.join(map(str, invalid))}")

    # Use provided probe_start / probe_end when available, otherwise fallback to pos ±100.
    fallback_start = (df["pos"] - 100).clip(lower=0)
    fallback_end = (df["pos"] + 100).clip(upper=df["chrom_size"])

    if "probe_start" in df.columns:
        df["probe_start"] = pd.to_numeric(df["probe_start"], errors="coerce")
        df["probe_start"] = df["probe_start"].fillna(fallback_start)
    else:
        df["probe_start"] = fallback_start

    if "probe_end" in df.columns:
        df["probe_end"] = pd.to_numeric(df["probe_end"], errors="coerce")
        df["probe_end"] = df["probe_end"].fillna(fallback_end)
    else:
        df["probe_end"] = fallback_end

    df[["probe_start", "probe_end"]] = df[["probe_start", "probe_end"]].astype(int)
    df["pos_id"] = df["chrom"].astype(str) + "_" + df["pos"].astype(str)
    df["chrom"] = pd.Categorical(
        df["chrom"], categories=chrom_df["chrom"].tolist(), ordered=True
    )
    df["pos_0"] = df["pos"] - 1
    return df.sort_values(by=["chrom", "pos"]).drop(columns=["chrom_size"])


def write_probe_id_file(
    df: pd.DataFrame, out_path: Path, probe_id: str
) -> pd.DataFrame:
    """
    写入 probe id 文件，并返回去重后的目标位点数据。
    """
    out_path.mkdir(parents=True, exist_ok=True)
    id_file = out_path / f"{probe_id}.id"
    unique_df = df.drop_duplicates(subset="pos_id")
    unique_df.to_csv(
        id_file,
        sep="\t",
        index=False,
        header=False,
        columns=["pos_id"],
    )
    return unique_df


def write_bed(df: pd.DataFrame, out_bed: Path, columns: Sequence[str]) -> None:
    """
    按 BED 三列格式写入文件。
    """
    df.to_csv(
        out_bed,
        sep="\t",
        index=False,
        header=False,
        columns=columns,
    )


def load_split_bed(split_bed: Path) -> pd.DataFrame:
    """
    读取并校验 split.bed 文件。
    """
    if not split_bed.exists():
        raise FileNotFoundError(f"找不到 split.bed: {split_bed}")

    try:
        split_bed_df = pd.read_table(
            split_bed,
            header=None,
            names=["chrom", "split_start", "split_end", "new_chrom"],
            usecols=[0, 1, 2, 3],
        )
    except Exception as exc:  # pragma: no cover - pandas 内部异常类型会随版本变化
        raise ValueError(
            "split.bed 格式错误，必须包含 4 列: chrom, split_start, split_end, new_chrom"
        ) from exc

    split_bed_df["chrom"] = split_bed_df["chrom"].astype(str)
    split_bed_df["new_chrom"] = split_bed_df["new_chrom"].astype(str)
    split_bed_df["split_start"] = pd.to_numeric(
        split_bed_df["split_start"], errors="raise"
    ).astype(int)
    split_bed_df["split_end"] = pd.to_numeric(
        split_bed_df["split_end"], errors="raise"
    ).astype(int)
    return split_bed_df


def load_split_chrom_order(split_genome_fai: Path) -> list[str]:
    """
    读取 split 基因组染色体顺序（来自 split.genome.fa.fai）。
    """
    if not split_genome_fai.exists():
        raise FileNotFoundError(f"找不到 split genome fai: {split_genome_fai}")

    chrom_df = pd.read_table(
        split_genome_fai,
        header=None,
        names=["chrom"],
        usecols=[0],
    )
    chrom_order = chrom_df["chrom"].astype(str).tolist()
    if not chrom_order:
        raise ValueError(f"split genome fai 为空: {split_genome_fai}")
    return chrom_order


def split_bed_dataframe(
    bed_df: pd.DataFrame,
    split_bed_df: pd.DataFrame,
    start_col: str,
    end_col: str,
) -> pd.DataFrame:
    """
    参考 gtf/split_bed.py，将 BED 坐标映射到 split 基因组坐标。
    """
    source_df = bed_df.copy()
    source_df["chrom"] = source_df["chrom"].astype(str)
    source_df["start"] = pd.to_numeric(source_df[start_col], errors="raise").astype(int)
    source_df["end"] = pd.to_numeric(source_df[end_col], errors="raise").astype(int)

    merge_df = source_df[["chrom", "start", "end"]].merge(split_bed_df, on="chrom")
    merge_df = merge_df[
        (merge_df["start"] >= merge_df["split_start"])
        & (merge_df["start"] < merge_df["split_end"])
    ].copy()
    merge_df["new_start"] = (merge_df["start"] - merge_df["split_start"]).astype(int)
    merge_df["new_end"] = (merge_df["end"] - merge_df["split_start"]).astype(int)
    return merge_df[["new_chrom", "new_start", "new_end"]]


def sort_split_output_by_fai(
    split_df: pd.DataFrame,
    split_chrom_order: list[str],
) -> pd.DataFrame:
    """
    按 split.genome.fa.fai 的染色体顺序排序 split 输出。
    """
    if split_df.empty:
        return split_df

    ordered_df = split_df.copy()
    ordered_df["new_chrom"] = ordered_df["new_chrom"].astype(str)
    missing = sorted(set(ordered_df["new_chrom"].unique()) - set(split_chrom_order))
    if missing:
        raise ValueError(
            f"split 输出中存在不在 split genome fai 的染色体: {', '.join(missing)}"
        )

    ordered_df["new_chrom"] = pd.Categorical(
        ordered_df["new_chrom"],
        categories=split_chrom_order,
        ordered=True,
    )
    ordered_df = ordered_df.sort_values(
        by=["new_chrom", "new_start", "new_end"],
        kind="mergesort",
    ).reset_index(drop=True)
    ordered_df["new_chrom"] = ordered_df["new_chrom"].astype(str)
    return ordered_df


def build_flank_intervals(
    df: pd.DataFrame, chrom_df: pd.DataFrame, flank_size: int
) -> pd.DataFrame:
    """
    Build flanking intervals for each probe.
    """
    df_for_merge = df.copy()
    df_for_merge["chrom"] = df_for_merge["chrom"].astype(str)
    merged = df_for_merge.merge(chrom_df, how="left")
    merged["chrom"] = pd.Categorical(
        merged["chrom"], categories=chrom_df["chrom"].tolist(), ordered=True
    )
    if merged["chrom_size"].isna().any():
        invalid = merged.loc[merged["chrom_size"].isna(), "chrom"].unique()
        raise ValueError(f"缺少染色体大小: {', '.join(map(str, invalid))}")

    merged[["flank_start", "flank_end"]] = merged.apply(
        lambda row: get_flank_start_end(
            row["probe_start"], row["probe_end"], flank_size, row["chrom_size"]
        ),
        axis=1,
        result_type="expand",
    )

    result = (
        merged.groupby("chrom", group_keys=True, sort=False, observed=False)
        .apply(
            lambda grp: merge_intervals(grp, "flank_start", "flank_end"),
            include_groups=False,
        )
        .reset_index(level=0)
        .reset_index(drop=True)
    )
    result[["flank_start", "flank_end"]] = result[["flank_start", "flank_end"]].astype(
        int
    )
    return result


@app.command(help=MODULE_HELP)
def main(
    design_table: Annotated[Path, typer.Argument(help="设计表路径（.tsv/.xlsx/.xls）")],
    genome_fai: Annotated[Path, typer.Argument(help="基因组 FASTA 索引（.fai）")],
    probe_id: Annotated[str, typer.Argument(help="输出文件名前缀")],
    out_path: Annotated[Path, typer.Argument(help="输出目录")],
    flank_size: Annotated[
        int, typer.Option(help="flanking 区域目标长度（必须大于 0）")
    ] = 200,
    split_bed: Annotated[
        Path | None,
        typer.Option(help="split.bed 文件路径；提供后仅输出拆分后的 *.bed"),
    ] = None,
    split_genome_fai: Annotated[
        Path | None,
        typer.Option(
            help="split 基因组 fai（.fai），用于按染色体顺序排序 split 输出；"
            "未提供时会尝试使用 split.bed 同目录下的 split.genome.fa.fai"
        ),
    ] = None,
) -> None:
    """
    根据设计表生成 panel 目标位点 bed / id 文件，并扩展 flanking 区域。
    """
    if flank_size <= 0:
        raise typer.BadParameter(f"flank_size 必须大于 0，当前值: {flank_size}")
    out_path.mkdir(parents=True, exist_ok=True)
    chrom_df = load_chrom_sizes(genome_fai)
    design_df = load_design_table(design_table)
    prepared_df = prepare_probe_dataframe(design_df, chrom_df)
    unique_df = write_probe_id_file(prepared_df, out_path, probe_id)
    flank_df = build_flank_intervals(prepared_df, chrom_df, flank_size)

    if split_bed is None:
        write_bed(
            unique_df,
            out_path / f"{probe_id}.bed",
            columns=["chrom", "pos_0", "pos"],
        )
        write_bed(
            flank_df,
            out_path / f"{probe_id}.snpcalling.bed",
            columns=["chrom", "flank_start", "flank_end"],
        )
        return

    split_bed_df = load_split_bed(split_bed)
    resolved_split_genome_fai = split_genome_fai
    if resolved_split_genome_fai is None:
        inferred_fai = split_bed.parent / "split.genome.fa.fai"
        if inferred_fai.exists():
            resolved_split_genome_fai = inferred_fai
        else:
            raise ValueError(
                "启用 --split-bed 时需要 --split-genome-fai，"
                "或在 split.bed 同目录提供 split.genome.fa.fai"
            )
    split_chrom_order = load_split_chrom_order(resolved_split_genome_fai)

    split_target_df = split_bed_dataframe(
        unique_df,
        split_bed_df,
        start_col="pos_0",
        end_col="pos",
    )
    split_target_df = sort_split_output_by_fai(split_target_df, split_chrom_order)
    split_snpcalling_df = split_bed_dataframe(
        flank_df,
        split_bed_df,
        start_col="flank_start",
        end_col="flank_end",
    )
    split_snpcalling_df = sort_split_output_by_fai(
        split_snpcalling_df,
        split_chrom_order,
    )
    write_bed(
        split_target_df,
        out_path / f"{probe_id}.bed",
        columns=["new_chrom", "new_start", "new_end"],
    )
    write_bed(
        split_snpcalling_df,
        out_path / f"{probe_id}.snpcalling.bed",
        columns=["new_chrom", "new_start", "new_end"],
    )


if __name__ == "__main__":
    app()
