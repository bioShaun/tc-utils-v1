import typer
import pandas as pd
from pathlib import Path


def get_score(row) -> float:
    score = 0
    if row["GC"] >= 0.35 and row["GC"] <= 0.65:
        score += 100
    allele_count = len(str(row["allele"]).split(","))
    isIndex = len(str(row["allele"])) > 3
    if isIndex:
        score -= 5
    if allele_count < 2:
        allele_count = 2
    if "maf" in row:
        if row["maf"] < 0.1:
            score -= 100
        if row["maf"] < 0.2:
            score -= 50
    return (
        score
        - 10 * (row["match_0"] - 1)
        + 1000 * (4 - row["score"])
        - (allele_count - 2) * 0.1
    )


def select_by_region(
    candidate_file: Path,
    region_cover: Path,
    region_count: int,
    output: Path,
    high_score: int = 0,
    label: str = "candidates",
) -> None:
    reg_df = pd.read_csv(
        region_cover, sep="\t", header=None, names=["region", "chrom", "pos"]
    )
    reg_df = reg_df[reg_df.chrom != "."]
    can_df = pd.read_csv(candidate_file)
    merged_df = reg_df.merge(can_df)
    merged_df["score2"] = merged_df.apply(get_score, axis=1)
    merged_df.sort_values("score2", ascending=False, inplace=True)
    center_df = merged_df[merged_df.sequence_type == "center"].copy()
    region_dfs = []
    for name, df in center_df.groupby("region"):
        high_score_df = df[df["score"] < high_score]
        low_score_df = df[df["score"] >= high_score]
        if len(high_score_df) >= region_count:
            region_dfs.append(high_score_df)
        else:
            region_dfs.append(high_score_df)
            region_dfs.append(low_score_df[: region_count - len(high_score_df)])
    final_df = pd.concat(region_dfs)
    final_df.sort_values(["chrom", "pos"], inplace=True)
    final_df["label"] = label
    final_df.to_csv(output, index=False, float_format="%.3f")


if __name__ == "__main__":
    typer.run(select_by_region)
