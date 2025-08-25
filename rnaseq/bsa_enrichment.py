from functools import partial
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import typer
from loguru import logger
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

GO_OUT_COLUMNS = [
    "Description",
    "Category",
    "GeneRatio",
    "BgRatio",
    "pvalue",
    "p.adjust",
    "geneID",
]

KEGG_OUT_COLUMNS = [
    "Description",
    "GeneRatio",
    "BgRatio",
    "pvalue",
    "p.adjust",
    "geneID",
]


def raw_hyper_test(total_gene: int, total_bg: int, row: pd.Series) -> float:
    k = row["pickCount"] - 1
    N = row["bgCount"]
    return hypergeom.sf(k, total_bg, total_gene, N)


def format_kegg_genes(gene_map_df: pd.DataFrame, gene_list: List[str]) -> str:
    out_list = []
    for gene in gene_list:
        ensemble_genes = gene_map_df.loc[gene]["ensembl"]
        if gene_map_df.loc[gene].shape[0] > 1:
            out_list.append(f"({','.join(ensemble_genes)})")
        else:
            out_list.append(ensemble_genes)
    return "/".join(out_list)


def format_df(enrich_df: pd.DataFrame, name_map_df: pd.DataFrame) -> pd.DataFrame:
    go2name_df = name_map_df.rename(
        columns={"name": "Description", "category": "Category"}
    )
    out_df = enrich_df.merge(go2name_df, left_index=True, right_index=True)
    out_df.index.name = "ID"
    out_df.sort_values(by="p.adjust", ascending=True, inplace=True)
    return out_df


def plot(enrich_df: pd.DataFrame, out_prefix: str) -> None:
    plot_data = enrich_df[:30].copy()
    plot_data["logPvalue"] = -np.log10(plot_data["p.adjust"])  # type: ignore
    plot_data = plot_data[plot_data["logPvalue"] > 0]
    plot_height = 2 + len(plot_data) * 0.2
    text_length = plot_data.Description.map(lambda x: len(x)).max()
    plot_width = 6 + text_length * 0.1
    sns.set(rc={"figure.figsize": (plot_width, plot_height)})
    sns.set_style("white")
    if "Category" in plot_data.columns:
        p = sns.barplot(
            y="Description", x="logPvalue", hue="Category", data=plot_data, dodge=False
        )
    else:
        p = sns.barplot(
            y="Description", x="logPvalue", data=plot_data, dodge=False, color="b"
        )
    p.set(xlabel="-log10(p.adjust)", ylabel="")
    plt.tight_layout()
    plt.savefig(f"{out_prefix}.png", dpi=300)
    plt.savefig(f"{out_prefix}.pdf")


def enrich_core(gene2go_df: pd.DataFrame, gene_list_df: pd.DataFrame):
    go2gene_df = pd.DataFrame(gene2go_df.groupby(["term"])["gene"].nunique()).rename(
        columns={"gene": "bgCount"}
    )
    current_gene2go_df = gene2go_df.merge(gene_list_df)
    current_go2gene_df = pd.DataFrame(
        current_gene2go_df.groupby(["term"])["gene"].nunique()
    ).rename(columns={"gene": "pickCount"})
    current_go2geneList_df = pd.DataFrame(
        current_gene2go_df.groupby(["term"])["gene"].unique()
    ).rename(columns={"gene": "geneID"})
    pickTotal = current_gene2go_df["gene"].nunique()
    bgTotal = gene2go_df["gene"].nunique()

    merged_term2gene_df = current_go2gene_df.merge(
        go2gene_df, left_index=True, right_index=True
    ).merge(current_go2geneList_df, left_index=True, right_index=True)
    hyper_test = partial(raw_hyper_test, pickTotal, bgTotal)
    merged_term2gene_df["pvalue"] = merged_term2gene_df.apply(hyper_test, axis=1)
    merged_term2gene_df["GeneRatio"] = [
        f"{each}|{pickTotal}" for each in merged_term2gene_df["pickCount"]
    ]
    merged_term2gene_df["BgRatio"] = [
        f"{each}|{bgTotal}" for each in merged_term2gene_df["bgCount"]
    ]
    merged_term2gene_df["p.adjust"] = multipletests(
        merged_term2gene_df["pvalue"], method="fdr_bh"
    )[1]
    return merged_term2gene_df


def go(gene_df: pd.DataFrame, id_map: Path, name_map: Path, out_prefix: str):
    gene2go_df = pd.read_csv(id_map)
    gene2go_df["gene"] = gene2go_df["gene"].astype(str)
    gene_list_df = gene_df.rename(columns={"Gene": "gene"}).drop_duplicates()
    gene_list_df["gene"] = gene_list_df["gene"].astype(str)
    enrich_df = enrich_core(gene2go_df, gene_list_df)
    name_map_df = pd.read_csv(name_map, sep="\t", index_col=0)
    out_df = format_df(enrich_df, name_map_df)
    out_df["geneID"] = out_df["geneID"].apply(lambda x: "/".join(x))
    out_df.to_csv(f"{out_prefix}.csv", columns=GO_OUT_COLUMNS)
    # plot(out_df, out_prefix)


def kegg(
    gene_df: pd.DataFrame,
    id_map: Path,
    name_map: Path,
    out_prefix: str,
    abbr: str,
    ncbi_map: Path = typer.Option(...),
):
    ncbi_map_df = pd.read_csv(ncbi_map, dtype={"ncbi": "str"})
    gene2kegg_df = pd.read_csv(id_map, dtype={"gene": "str"})
    ens_gene_list_df = gene_df.rename(columns={"Gene": "gene"}).drop_duplicates()
    ncbi_gene_df = ncbi_map_df.merge(
        ens_gene_list_df, left_on="ensembl", right_on="gene"
    )
    ncbi_gene_list_df = (
        ncbi_gene_df[["ncbi"]].drop_duplicates().rename(columns={"ncbi": "gene"})
    )
    enrich_df = enrich_core(gene2kegg_df, ncbi_gene_list_df)
    name_map_df = pd.read_csv(name_map, dtype={"term": "str"})
    name_map_df["term"] = [f"{abbr}{each}" for each in name_map_df["term"]]
    name_map_df = name_map_df.set_index("term")
    out_df = format_df(enrich_df, name_map_df)
    ncbi_map_df = ncbi_map_df.set_index("ncbi")
    current_format_id = partial(format_kegg_genes, ncbi_map_df)
    out_df["geneID"] = out_df["geneID"].map(current_format_id)
    out_df.to_csv(f"{out_prefix}.csv", columns=KEGG_OUT_COLUMNS)
    # plot(out_df, out_prefix)


def bsa_enrichment(bsa_dir: Path, ann_dir: Path, abbr: str) -> None:
    kegg_id_map = ann_dir / "kegg.idmap.csv"
    kegg_name = ann_dir / "kegg.name.csv"
    go_name = ann_dir / "go.name.csv"
    go_id_map = ann_dir / "go.idmap.csv"
    ncbi_id_map = ann_dir / "ncbi.idmap.csv"
    for csv_file in bsa_dir.rglob("*.csv"):
        if "%.csv" in csv_file.name or "CI" in csv_file.name:
            if csv_file.parent.name == "gene_enrichment":
                continue
            df = pd.read_csv(csv_file)
            gene_df = df[["Gene"]].copy()
            out_dir = csv_file.parent / "gene_enrichment"
            out_dir.mkdir(parents=True, exist_ok=True)
            go_out_prefix = out_dir / f"{csv_file.stem}.go"
            kegg_out_prefix = out_dir / f"{csv_file.stem}.kegg"
            logger.info(f"Running GO enrichment for {csv_file}")
            try:
                go(
                    gene_df=gene_df,
                    id_map=go_id_map,
                    name_map=go_name,
                    out_prefix=go_out_prefix,
                )
            except Exception as e:
                logger.error(f"Error running GO enrichment for {csv_file}")
                logger.error(e)
            logger.info(f"Running KEGG enrichment for {csv_file}")
            try:
                kegg(
                    gene_df=gene_df,
                    id_map=kegg_id_map,
                    name_map=kegg_name,
                    out_prefix=kegg_out_prefix,
                    abbr=abbr,
                    ncbi_map=ncbi_id_map,
                )
            except Exception as e:
                logger.error(f"Error running KEGG enrichment for {csv_file}")
                logger.error(e)


if __name__ == "__main__":
    typer.run(bsa_enrichment)
