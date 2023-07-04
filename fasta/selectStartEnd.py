import typer
import pandas as pd

from pathlib import Path


def selectStartEnd(seq_tab: Path, out_file: Path) -> None:
    seq_df = pd.read_csv(
        seq_tab,
        sep="\t",
        header=None,
        names=["slide_name", "sequence", "empty_col", "seq_length", "gc_content"],
    )

    seq_df["name"] = seq_df["slide_name"].map(lambda x: x.split("_")[0])
    start_seq_df = seq_df.groupby(["name"]).first()
    end_seq_df = seq_df.groupby(["name"]).last()
    start_end_seq_df = pd.concat([start_seq_df, end_seq_df])
    start_end_seq_df = start_end_seq_df.reset_index()
    start_end_seq_df.to_csv(
        out_file,
        sep="\t",
        header=False,
        columns=["slide_name", "sequence", "gc_content"],
        index=False,
    )


if __name__ == "__main__":
    typer.run(selectStartEnd)
