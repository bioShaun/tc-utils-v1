from typing import List, Optional
import typer

import pandas as pd

from pathlib import Path
from typing_extensions import Annotated


def main(
    candidates: Path,
    to_remove: Path,
    new_candidates: Path,
    to_add: Optional[List[Path]] = typer.Option(None),
) -> None:
    can_df = pd.read_csv(candidates)
    to_remove_df = pd.read_csv(
        to_remove, header=None, names=["chrom", "pos", "allele", "score"]
    )
    to_remove_df["id"] = to_remove_df.apply(
        lambda x: f"{x['chrom']}_{x['pos']}_center", axis=1
    )
    remove_can_df = can_df[~can_df["id"].isin(to_remove_df["id"])].copy()
    if to_add is None:
        new_can_df = remove_can_df
        new_can_df.sort_values(["chrom", "pos"], inplace=True)
        new_can_df.to_csv(new_candidates, index=False)
    else:
        to_add_df = pd.concat([pd.read_csv(each) for each in to_add])
        # to_add_df = pd.read_csv(to_add)
        new_can_df = pd.concat([remove_can_df, to_add_df])
        new_can_df.sort_values(["chrom", "pos"], inplace=True)
        new_can_df.to_csv(new_candidates, index=False)


if __name__ == "__main__":
    typer.run(main)
