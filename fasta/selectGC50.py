from pathlib import Path
import typer
import pandas as pd


def chooseGC50(seq_tab: Path, out_file: Path) -> None:
    seq_df = pd.read_csv(
        seq_tab,
        sep="\t",
        header=None,
        names=["slide_name", "sequence", "empty_col", "seq_length", "gc_content"],
    )

    seq_df["name"] = seq_df["slide_name"].map(lambda x: x.split("_")[0])
    seq_df["gc_50_abs"] = seq_df["gc_content"].map(lambda x: abs(x - 50))
    best_seq_df = seq_df.take(seq_df.groupby(["name"])["gc_50_abs"].idxmin())
    best_seq_df.to_csv(
        out_file,
        sep="\t",
        header=False,
        columns=["name", "sequence", "gc_content"],
        index=False,
    )


if __name__ == "__main__":
    typer.run(chooseGC50)
