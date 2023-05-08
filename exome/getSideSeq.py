from pathlib import Path
import typer
import pandas as pd


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
    return (
        score
        - 10 * (row["match_0"] - 1)
        + 1000 * (4 - row["score"])
        - (allele_count - 2) * 0.1
    )


def candidates_side_seq(
    selected_file: Path, candidate_file: Path, output: Path, report_best: bool = True
) -> None:
    select_df = pd.read_csv(selected_file)
    select_loci = select_df[["chrom", "pos"]].drop_duplicates()
    can_df = pd.read_csv(candidate_file)
    all_side_df = can_df[can_df["sequence_type"] != "center"]
    select_side_seqs = all_side_df.merge(select_loci)
    select_side_seqs["score2"] = select_side_seqs.apply(get_score, axis=1)
    if report_best:
        select_side_seqs.sort_values("score2", ascending=False, inplace=True)
        select_side_seqs.drop_duplicates(["chrom", "pos"], inplace=True)
    select_side_seqs.sort_values(["chrom", "pos"], inplace=True)
    select_side_seqs.to_csv(output, index=False)


if __name__ == "__main__":
    typer.run(candidates_side_seq)
