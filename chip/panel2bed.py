"""
设计表生成 panel 目标位点 bed / id 文件，并扩展 flanking 区域。
"""

from pathlib import Path
from typing import Iterable, Tuple

import pandas as pd
import typer

REQUIRED_COLUMNS = ("chrom", "pos")


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


def write_probe_targets(df: pd.DataFrame, out_path: Path, probe_id: str) -> None:
    """
    Write probe targets to a bed file and an id file.
    """
    out_path.mkdir(parents=True, exist_ok=True)
    id_file = out_path / f"{probe_id}.id"
    target_bed_file = out_path / f"{probe_id}.bed"
    unique_df = df.drop_duplicates(subset="pos_id")
    unique_df.to_csv(
        id_file,
        sep="\t",
        index=False,
        header=False,
        columns=["pos_id"],
    )
    unique_df.to_csv(
        target_bed_file,
        sep="\t",
        index=False,
        header=False,
        columns=["chrom", "pos_0", "pos"],
    )


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


def main(
    design_table: Path,
    genome_fai: Path,
    probe_id: str,
    out_path: Path,
    flank_size: int = 200,
) -> None:
    """
    根据设计表生成 panel 目标位点 bed / id 文件，并扩展 flanking 区域。

    参数:
        design_table: 设计表（tsv/xlsx）
        genome_fai: genome fasta index，用于获取染色体大小
        probe_id: 输出文件名前缀
        out_path: 输出目录
        flank_size: flanking 区域目标长度
    """
    chrom_df = load_chrom_sizes(genome_fai)
    design_df = load_design_table(design_table)
    prepared_df = prepare_probe_dataframe(design_df, chrom_df)
    write_probe_targets(prepared_df, out_path, probe_id)
    flank_df = build_flank_intervals(prepared_df, chrom_df, flank_size)
    flank_df.to_csv(
        out_path / f"{probe_id}.snpcalling.bed", sep="\t", index=False, header=False
    )


if __name__ == "__main__":
    typer.run(main)
