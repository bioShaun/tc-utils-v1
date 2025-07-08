from pathlib import Path

import polars as pl
import typer
from loguru import logger


def fetch_lazy_df(data_dir: Path, pattern: str) -> pl.LazyFrame:
    return pl.concat([pl.scan_csv(each) for each in data_dir.glob(pattern)])


def qtlseqr_snpIndex_filter(df: pl.LazyFrame, ci_col: str) -> pl.LazyFrame:
    snpindex_col = "tricubeDeltaSNP"
    return df.filter(
        (pl.col(snpindex_col) > abs(pl.col(ci_col)))
        | (pl.col(snpindex_col) < -abs(pl.col(ci_col)))
    )


def qtlseqr_candidate_filter(df: pl.LazyFrame, stats_col: str, ci: str) -> pl.LazyFrame:
    df_stats = (
        df.select(
            [
                pl.col(stats_col).mean().alias("mean"),
                pl.col(stats_col).std().alias("std"),
            ]
        )
        .collect()
        .to_dict(as_series=False)
    )

    mean_stats = df_stats["mean"][0]
    std_stats = df_stats["std"][0]

    logger.info(f"mean_stats: {mean_stats}, std_stats: {std_stats}")

    if ci == "CI_95":
        ci_factor = 1.96
    elif ci == "CI_99":
        ci_factor = 2.58
    else:
        raise ValueError(f"Invalid CI: {ci}")

    ci_upper = mean_stats + ci_factor * std_stats
    logger.info(f"ci_upper: {ci_upper}")

    df_filtered = df.filter((pl.col(stats_col) >= ci_upper))
    return df_filtered


def varbscore_candidate_filter(data_dir: Path, ci: str) -> pl.DataFrame:
    df = fetch_lazy_df(data_dir, "*varBScore/data/*.csv")
    stats_col = "varBScore"
    varbscore_unique = df.unique(subset=["CHROM", "Start", "End"]).collect()

    stats_df = varbscore_unique.select(
        [
            pl.col("varBScore").mean().alias("mean"),
            pl.col("varBScore").std().alias("std"),
        ]
    ).to_dict(as_series=False)

    mean_ed = stats_df["mean"][0]
    std_ed = stats_df["std"][0]

    if ci == "CI_95":
        ci_factor = 1.96
    elif ci == "CI_99":
        ci_factor = 2.58
    else:
        raise ValueError(f"Invalid CI: {ci}")

    ci_upper = mean_ed + ci_factor * std_ed

    df_filtered = df.filter(pl.col(stats_col) >= ci_upper)
    return df_filtered.collect().sort(["CHROM", "Start", "End"])


def qtlseqr_pos_label(
    gprime_df: pl.LazyFrame, ed_df: pl.LazyFrame, snpindex_df: pl.LazyFrame
) -> pl.DataFrame:
    ed_df_ci95_pos = ed_df.select(["CHROM", "POS"]).with_columns(
        pl.lit("ED").alias("Label")
    )
    gprime_df_ci95_pos = gprime_df.select(["CHROM", "POS"]).with_columns(
        pl.lit("Gprime").alias("Label")
    )
    snpindex_df_ci95_pos = snpindex_df.select(["CHROM", "POS"]).with_columns(
        pl.lit("SNPindex").alias("Label")
    )
    all_pos_df = pl.concat(
        [ed_df_ci95_pos, gprime_df_ci95_pos, snpindex_df_ci95_pos]
    ).collect()
    return (
        all_pos_df.group_by(["CHROM", "POS"])
        .agg(pl.col("Label").sort().unique().sort().alias("Label"))
        .with_columns(pl.col("Label").list.join(", "))
    )


def merge_qtlseqr_candidate(data_dir: Path, ci: str) -> pl.DataFrame:
    qtlseqr_df = fetch_lazy_df(data_dir, "*QTLseqr*/data/*.csv")
    snpindex_df = qtlseqr_snpIndex_filter(qtlseqr_df, ci)
    gprime_df = qtlseqr_candidate_filter(qtlseqr_df, "Gprime", ci)
    ed_df = qtlseqr_candidate_filter(qtlseqr_df, "fitted", ci)
    merged_df = (
        pl.concat([snpindex_df, gprime_df, ed_df])
        .unique(subset=["CHROM", "POS"])
        .collect()
    )

    pos_label_df = qtlseqr_pos_label(gprime_df, ed_df, snpindex_df)
    return merged_df.join(pos_label_df, on=["CHROM", "POS"], how="left").sort(
        ["CHROM", "POS"]
    )


def fetch_qtlseqr_dir(data_dir: Path) -> Path:
    for each_path in data_dir.glob("*QTLseqr*"):
        if each_path.is_dir():
            return each_path
    raise FileNotFoundError(f"QTLseqr dir not found in {data_dir}")


def fetch_varbscore_dir(data_dir: Path) -> Path:
    for each_path in data_dir.glob("*varBScore*"):
        if each_path.is_dir():
            return each_path
    raise FileNotFoundError(f"VarBscore dir not found in {data_dir}")


def main(
    bsa_dir: Path = typer.Argument(..., help="BSA dir"),
) -> None:
    for ci in ["CI_95", "CI_99"]:
        # VarBscore candidate
        varbscore_ci_df = varbscore_candidate_filter(bsa_dir, ci)
        varbscore_dir = fetch_varbscore_dir(bsa_dir)
        varbscore_ci_file = varbscore_dir / f"VarBscore_{ci}.csv"
        varbscore_ci_df.write_csv(varbscore_ci_file)

        # QTLseqr candidate
        qtlseqr_ci_df = merge_qtlseqr_candidate(bsa_dir, ci)
        qtlseqr_dir = fetch_qtlseqr_dir(bsa_dir)
        qtlseqr_ci_file = qtlseqr_dir / f"QTLseqr_{ci}.csv"
        qtlseqr_ci_df.write_csv(qtlseqr_ci_file)


if __name__ == "__main__":
    typer.run(main)
