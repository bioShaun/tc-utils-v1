from pathlib import Path

import pandas as pd
import typer

FLANK_SIZE = 60


def main(paf: Path, id_map: Path) -> None:
    paf_df = pd.read_table(paf, header=None)
    paf_df.sort_values([11, 9], ascending=False, inplace=True)
    paf_df.drop_duplicates(subset=[0], inplace=True)
    filter_df = paf_df[paf_df[9] > FLANK_SIZE].copy()
    filter_df["pos"] = filter_df[7] + 60 + 1 - filter_df[2]
    filter_df["id"] = filter_df.apply(lambda x: f'{x[5]}_{x["pos"]}', axis=1)
    filter_df.rename(columns={0: "old_id"}, inplace=True)
    filter_df.to_csv(
        id_map, header=False, index=False, columns=["old_id", "id"], sep="\t"
    )


if __name__ == "__main__":
    typer.run(main)
