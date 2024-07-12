from pathlib import Path

import pandas as pd
import typer


def main(mean_q: Path, sample_list: Path, out_file: Path) -> None:
    sample_df = pd.read_csv(sample_list, header=None, names=["sample_id"])
    mean_q_df = pd.read_table(mean_q, sep="\s+", header=None)
    mean_q_df.columns = [f"pop{each+1}" for each in mean_q_df.columns]
    merged_df = pd.concat([sample_df, mean_q_df], axis=1)
    melt_df = merged_df.melt(id_vars="sample_id", var_name="group")
    sample_group_df = melt_df.loc[melt_df.groupby(["sample_id"])["value"].idxmax()]
    sample_group_df.to_csv(
        out_file, sep="\t", index=False, columns=["sample_id", "group"]
    )


if __name__ == "__main__":
    typer.run(main)
