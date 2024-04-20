from pathlib import Path

import pandas as pd
import typer

FLANK_SIZE = 60


def main(paf: Path, id_map: Path) -> None:
    paf_df = pd.read_table(
        paf,
        header=None,
        usecols=[0, 2, 5, 7, 9, 11],
        names=["id", "probe_start", "chrom", "match_start", "match_length", "mapq"],
    )
    paf_df.sort_values(["mapq", "match_length"], ascending=False, inplace=True)
    paf_df.drop_duplicates(subset=["id"], inplace=True)
    filter_df = paf_df[paf_df["match_length"] > FLANK_SIZE].copy()
    filter_df["pos"] = filter_df["match_start"] + 60 + 1 - filter_df["probe_start"]
    filter_df["new_id"] = filter_df.apply(lambda x: f'{x["chrom"]}_{x["pos"]}', axis=1)
    filter_df.to_csv(
        id_map, header=False, index=False, columns=["id", "new_id"], sep="\t"
    )


if __name__ == "__main__":
    typer.run(main)
