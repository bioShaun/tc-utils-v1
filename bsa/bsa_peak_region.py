from pathlib import Path
from typing import Optional

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


def varfilter_candidate_regions(data_dir: Path, ci_number: float) -> pl.DataFrame:
    df = fetch_lazy_df(data_dir, "*varFilter/var.filter.density.csv")
    df_stats = df.select(
        [pl.col("variantCount").quantile(ci_number).alias("q_value")]
    ).collect()

    q_value = df_stats["q_value"][0]
    df_filtered = df.filter(pl.col("variantCount") > q_value).collect()
    return df_filtered


def varfilter_candidate_filter(data_dir: Path, ci_number: float) -> pl.DataFrame:
    candidate_region = varfilter_candidate_regions(data_dir, ci_number)
    candidate_df_list = []
    varfilter_df = fetch_lazy_df(data_dir, "*varFilter/data/*.csv")
    for row in candidate_region.iter_rows(named=True):
        candidate_df_list.append(
            varfilter_df.filter(
                (pl.col("CHROM") == row["chrom"])
                & (pl.col("POS") >= row["start"])
                & (pl.col("POS") <= row["end"])
            )
        )
    return pl.concat(candidate_df_list).collect()


def qtlseqr_candidate_filter(df: pl.LazyFrame, stats_col: str, ci: str) -> pl.LazyFrame:
    df_unique = df.unique(subset=["CHROM", "POS"])
    df_stats = (
        df_unique.select(
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
        all_pos_df.groupby(["CHROM", "POS"])
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
        .unique(subset=["CHROM", "POS", "REF", "ALT", "Feature", "Transcript"])
        .collect()
    )

    pos_label_df = qtlseqr_pos_label(gprime_df, ed_df, snpindex_df)
    return merged_df.join(pos_label_df, on=["CHROM", "POS"], how="left").sort(
        ["CHROM", "POS"]
    )


def fetch_analysis_dir(data_dir: Path, pattern: str) -> Optional[Path]:
    for each_path in data_dir.glob(pattern):
        if each_path.is_dir():
            return each_path
    return None


def main(
    bsa_dir: Path = typer.Argument(..., help="BSA dir"),
) -> None:
    for ci in ["CI_95", "CI_99"]:
        # VarFilter candidate
        varFilter_dir = fetch_analysis_dir(bsa_dir, "*varFilter*")
        if varFilter_dir is not None:
            ci_number = int(ci.split("_")[1])
            q_number = ci_number / 100
            varfilter_ci_df = varfilter_candidate_filter(bsa_dir, q_number)
            varfilter_ci_file = varFilter_dir / f"VarFilter_top{100 - ci_number}%.csv"
            varfilter_ci_df.write_csv(varfilter_ci_file)
        # VarBscore candidate
        varbscore_dir = fetch_analysis_dir(bsa_dir, "*varBScore*")
        if varbscore_dir is not None:
            varbscore_ci_df = varbscore_candidate_filter(bsa_dir, ci)
            varbscore_ci_file = varbscore_dir / f"VarBscore_{ci}.csv"
            varbscore_ci_df.write_csv(varbscore_ci_file)

        # QTLseqr candidate
        qtlseqr_dir = fetch_analysis_dir(bsa_dir, "*QTLseqr*")
        if qtlseqr_dir is not None:
            qtlseqr_ci_df = merge_qtlseqr_candidate(bsa_dir, ci)
            qtlseqr_ci_file = qtlseqr_dir / f"QTLseqr_{ci}.csv"
            qtlseqr_ci_df.write_csv(qtlseqr_ci_file)


if __name__ == "__main__":
    typer.run(main)
