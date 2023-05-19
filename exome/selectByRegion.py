import typer
import pandas as pd
from pathlib import Path


def get_score(row) -> float:
    score = 0
    if row["GC"] >= 0.35 and row["GC"] <= 0.65:
        score += 100
    allele_count = len(str(row["allele"]).split(","))
    isIndel = len(str(row["allele"])) > 3
    if isIndel:
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
    output_prefix: Path,
    auto_add: bool = True,
    region_select_count: int = typer.Option(None),
    region_select_bed: Path = typer.Option(None),
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
    if region_select_bed:
        region_select_df = pd.read_csv(
            region_select_bed, sep="\t", header=None, names=["region", "count"]
        )
    elif region_select_count:
        region_select_df = reg_df[["region"]]
        region_select_df["count"] = region_select_count
        region_select_df.drop_duplicates(inplace=True)
    else:
        raise ValueError("region_select_count and region_select_bed both none")
    missed_count = 0
    missed_region_list = []
    filter_df_list = []
    for _, region_df in region_select_df.iterrows():
        df = center_df[center_df.region == region_df.region]
        high_score_df = df[df["score"] <= high_score]
        low_score_df = df[df["score"] > high_score]
        region_count = region_df["count"]
        if len(df) < region_count:
            missed_count += region_count - len(df)
            missed_region_list.append([region_df.region, region_count - len(df)])
        if len(df) > region_count:
            filter_df_list.append(df[region_count:])
        if len(high_score_df) >= region_count:
            region_dfs.append(high_score_df)
        else:
            region_dfs.append(high_score_df)
            region_dfs.append(low_score_df[: region_count - len(high_score_df)])
    final_df = pd.concat(region_dfs)
    print(missed_count)
    if auto_add:
        filter_out_df = center_df[~center_df["id"].isin(final_df["id"])].copy()
        print(len(filter_out_df))
        tmp_df = filter_out_df[:missed_count]
        minimal_score = tmp_df["score2"].min()
        add_df = filter_out_df[filter_out_df["score2"] >= minimal_score].copy()
        add_df = add_df.sample(n=missed_count, random_state=1)
        auto_add_file = output_prefix.with_suffix(".auto_add.csv")
        add_df.to_csv(auto_add_file, index=False, float_format="%.3f")
        final_df = pd.concat([final_df, add_df])
    final_df.sort_values(["chrom", "pos"], inplace=True)
    candidate_file = output_prefix.with_suffix(".select.csv")
    missed_region_df = pd.DataFrame(missed_region_list, columns=["region", "count"])
    missed_region_file = output_prefix.with_suffix(".missed_region.csv")
    missed_region_df.to_csv(missed_region_file, index=False, sep="\t", header=False)
    final_df["label"] = label
    final_df.to_csv(candidate_file, index=False, float_format="%.3f")


if __name__ == "__main__":
    typer.run(select_by_region)
