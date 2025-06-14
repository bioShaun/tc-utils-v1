from pathlib import Path
from typing import List, Tuple

import gtfparse
import pandas as pd
import typer


def gtf2gene_locations(gtf_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Parses a GTF file and extracts the gene location information.

    Args:
        gtf_path (Path): The path to the GTF file.

    Returns:
        pd.DataFrame: A DataFrame containing the gene location information.
    """
    gtf_df = gtfparse.read_gtf(gtf_path, result_type="pandas")
    if isinstance(gtf_df, pd.DataFrame):
        genes = gtf_df.groupby("gene_id")
        tr2gene = gtf_df[["transcript_id", "gene_id"]].dropna().drop_duplicates()
        locations = pd.DataFrame(
            {
                "Chrom": genes["seqname"].first(),
                "Start": genes["start"].min(),
                "End": genes["end"].max(),
                "Strand": genes["strand"].first(),
            }
        )
        # locations.loc[~locations["Strand"].isin(["+", "-"]), "Strand"] = "."
        locations = locations.astype(str)
        locations["Locus"] = (
            locations["Chrom"]
            + ":"
            + locations["Start"]
            + "-"
            + locations["End"]
            + "|"
            + locations["Strand"]
        )
        locations = locations.reset_index()[["gene_id", "Locus"]]
        return tr2gene, locations
    raise ValueError(f"{gtf_path} is not a valid GTF file")


def go_ko_df(gene_eggnog: pd.DataFrame, name: str) -> pd.DataFrame:
    term_df = gene_eggnog[["gene_id", name]].copy()
    term_df = term_df[term_df[name] != "-"]
    term_df[name] = term_df[name].map(lambda x: x.split(","))
    term_df = term_df.explode(name)
    term_df.drop_duplicates(inplace=True)
    term_df.columns = ["gene", "term"]
    return term_df


def uniq_gene_des(names: List[str]) -> str:
    if len(names) == 1:
        return names[0]
    non_na_names = [name for name in names if name != "-"]
    return "|".join(non_na_names)


def main(
    eggnog_file: Path,  # type: Path
    gtf: str,  # type: str
    output_dir: Path,  # type: Path
) -> None:
    """
    Annotates gene information from eggNOGmapper output.

    Args:
        eggnog_file: eggNOGmapper output file.
        gtf: GTF file used to obtain gene location information.
        output_dir: Directory where to save output files.

    Returns:
        None
    """
    eggnog_df = pd.read_table(eggnog_file, header=None, comment="#")
    tr2gene, gene_locations = gtf2gene_locations(gtf)
    eggnog_df = eggnog_df[[0, 8, 7, 9, 12]].copy()
    eggnog_df.columns = ["transcript_id", "gene_name", "gene_description", "go", "ko"]
    gene_egg = tr2gene.merge(eggnog_df, how="left")
    gene_egg.fillna("-", inplace=True)
    gene_egg.drop_duplicates(subset=["transcript_id"], inplace=True)
    gene_egg.drop("transcript_id", axis=1, inplace=True)

    non_uniq_gene_description = gene_egg[["gene_id", "gene_name", "gene_description"]]
    gene_names = (
        non_uniq_gene_description.groupby(["gene_id"])["gene_name"]
        .unique()
        .map(uniq_gene_des)
    )
    gene_des = (
        non_uniq_gene_description.groupby(["gene_id"])["gene_description"]
        .unique()
        .map(uniq_gene_des)
    )
    rm_dup_des_df = pd.concat([gene_names, gene_des], axis=1).reset_index()
    gene_description = gene_locations.merge(rm_dup_des_df)
    gene_description.to_csv(output_dir / "gene.description.txt", index=False, sep="\t")

    go_id_map = go_ko_df(gene_egg, "go")
    go_id_map.to_csv(output_dir / "go.idmap.csv", index=False)
    ko_id_map = go_ko_df(gene_egg, "ko")
    ko_id_map = ko_id_map[ko_id_map["term"].str.startswith("ko")]
    ko_id_map.to_csv(output_dir / "kegg.idmap.csv", index=False)
    kegg_gene_map = ko_id_map[["gene", "gene"]]
    kegg_gene_map.columns = ["ensembl", "ncbi"]
    kegg_gene_map.to_csv(output_dir / "ncbi.idmap.csv", index=False)


if __name__ == "__main__":
    typer.run(main)
