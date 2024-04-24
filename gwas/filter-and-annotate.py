import typer
import pandas as pd
from pathlib import Path

from tqdm import tqdm
from loguru import logger

TESTS = ["wald", "lrt", "score"]
PVALUE = [1e-5, 5e-8]
PVALUE_TEXT = ["1e-5", "5e-8"]


def filter_by_pvalue(
    df: pd.DataFrame,
    pvalue: float,
) -> pd.DataFrame:
    """Filters the dataframe by p-value."""
    masks = [df[f"p_{test}"] <= pvalue for test in TESTS]
    mask = pd.concat(masks, axis=1).any(axis=1)
    return df[mask].drop_duplicates(ignore_index=True)


def main(output_dir: Path, annotation_file: Path) -> None:
    """
    Filters GWAS results and annotates the remaining results.

    Reads association files from the output directory and merges them with
    variant information from the annotation file. The resulting dataframe is
    saved to a file with the suffix ".anno.txt". The function also filters the
    results by p-value and saves the filtered results to files with the suffix
    ".PVALUE_TEXT.anno.txt", where PVALUE_TEXT is a string representing the
    p-value threshold.

    Parameters:
        output_dir (Path): The directory containing the association files.
        annotation_file (Path): The path to the variant annotation file.

    Returns:
        None
    """
    # Create a dataframe with variant information.
    annotation_df = pd.read_table(annotation_file)
    annotation_df["rs"] = annotation_df.apply(
        lambda row: f"{row['chrom']}_{row['pos']}", axis=1
    )
    annotation_df.drop(columns=["chrom", "pos"], inplace=True)

    # Iterate over the association files in the output directory.
    for association_file in tqdm(output_dir.glob("*.assoc.txt")):
        # Read the association file and merge it with the annotation file.
        logger.info(f"annotating: {association_file}...")
        df = pd.read_table(association_file)
        df = df.merge(annotation_df, on="rs")
        df["transcript_start"] = df["transcript_start"].astype("Int64")
        df["transcript_end"] = df["transcript_end"].astype("Int64")

        # Save the annotated results to a file.
        annotated_file = association_file.with_suffix(".anno.txt")
        df.to_csv(annotated_file, sep="\t", index=False)

        # Filter the results by p-value and save the filtered results.
        for n, pvalue in enumerate(PVALUE):
            filtered_df = filter_by_pvalue(df, pvalue)
            filtered_file = association_file.with_suffix(f".{PVALUE_TEXT[n]}.anno.txt")
            filtered_df.to_csv(filtered_file, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)
