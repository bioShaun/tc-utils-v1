import tempfile
from pathlib import Path

import pandas as pd
import typer

LOCATION_COLS = ["CHROM", "POS", "REF", "ALT"]

COLUMN_MAP = {
    "sample": "样本",
    "total_count": "群体位点数",
    "missing_count": "缺失位点数",
    "missing_fraction": "缺失率",
    "fraction": "检出率",
}


def test_main():
    # Create temporary files
    dp_file = tempfile.NamedTemporaryFile(delete=False)
    sample_file = tempfile.NamedTemporaryFile(delete=False)
    out_file_prefix = tempfile.NamedTemporaryFile(delete=False).name

    # Write test data to temporary files
    dp_data = """1	100	A	T	1	0	.
        1	200	G	C	0	1	1
        1	300	T	A	1	1	0"""
    sample_data = """sample1
sample2
sample3"""
    dp_file.write(dp_data.encode())
    dp_file.close()
    sample_file.write(sample_data.encode())
    sample_file.close()

    # Call the main function
    main(Path(dp_file.name), Path(sample_file.name), Path(out_file_prefix))

    # Read the output files and check their contents
    depth_df = pd.read_excel(f"{out_file_prefix}_depth.xlsx")
    stats_df = pd.read_excel(f"{out_file_prefix}_stats.xlsx")
    depth_df["CHROM"] = depth_df["CHROM"].astype("str")
    expected_depth_df = pd.DataFrame(
        {
            "CHROM": ["1", "1", "1"],
            "POS": [100, 200, 300],
            "REF": ["A", "G", "T"],
            "ALT": ["T", "C", "A"],
            "sample1": [1, 0, 1],
            "sample2": [0, 1, 1],
            "sample3": [0, 1, 0],
        }
    )
    expected_stats_df = pd.DataFrame(
        {
            "sample": ["sample1", "sample2", "sample3"],
            "total_count": [3, 3, 3],
            "missing_count": [1, 1, 2],
            "missing_fraction": [1 / 3, 1 / 3, 2 / 3],
            "fraction": [2 / 3, 2 / 3, 1 / 3],
        }
    ).rename(columns=COLUMN_MAP)

    pd.testing.assert_frame_equal(depth_df, expected_depth_df)
    pd.testing.assert_frame_equal(stats_df, expected_stats_df)

    # Clean up temporary files
    import os

    os.remove(dp_file.name)
    os.remove(sample_file.name)
    os.remove(f"{out_file_prefix}_depth.xlsx")
    os.remove(f"{out_file_prefix}_stats.xlsx")


def main(dp_file: Path, sample_file: Path, out_file_prefix: Path) -> None:
    samples = pd.read_table(sample_file, header=None)[0].tolist()
    print(samples)
    dp_df = pd.read_table(
        dp_file, header=None, names=[*LOCATION_COLS, *samples]
    ).set_index(LOCATION_COLS)
    print(dp_df)
    dp_df.replace(".", 0, inplace=True)
    dp_df = dp_df.astype("int")
    passed_df = dp_df.gt(0).sum().reset_index()
    passed_df.columns = ["sample", "count"]
    passed_df["fraction"] = passed_df["count"] / dp_df.shape[0]
    passed_df["missing_fraction"] = 1 - passed_df["fraction"]
    passed_df["total_count"] = dp_df.shape[0]
    passed_df["missing_count"] = passed_df["total_count"] - passed_df["count"]
    dp_df = dp_df.reset_index()
    dp_df.to_excel(f"{out_file_prefix}_depth.xlsx", index=False)
    passed_df.rename(columns=COLUMN_MAP, inplace=True)
    passed_df.to_excel(
        f"{out_file_prefix}_stats.xlsx",
        index=False,
        columns=[
            "样本",
            "群体位点数",
            "缺失位点数",
            "缺失率",
            "检出率",
        ],
    )


if __name__ == "__main__":
    typer.run(main)
