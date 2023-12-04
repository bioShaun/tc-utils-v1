from pathlib import Path
import typer
import pandas as pd


def parse_kegg(kegg_term: str) -> str:
    kegg_id = kegg_term.split(": ")[-1]
    return f'ko{kegg_id.split("+")[0]}'


def main(interpro_file: Path, gene2tr: Path, out_dir: Path) -> None:
    interpro_info_list = []
    with interpro_file.open() as interpro_inf:
        for eachline in interpro_inf:
            eachline_list = eachline.strip().split("\t")
            go = None
            kegg = None
            tr_id = eachline_list[0]
            if len(eachline_list) >= 14:
                if eachline_list[13]:
                    go = eachline_list[13].split("|")
            if len(eachline_list) >= 15:
                if eachline_list[14]:
                    kegg_list = [
                        parse_kegg(each)
                        for each in eachline_list[14].split("|")
                        if each.startswith("KEGG")
                    ]
                    if kegg_list:
                        kegg = kegg_list
            interpro_info_list.append({"tr_id": tr_id, "go": go, "kegg": kegg})
    interpro_df = pd.DataFrame(interpro_info_list)
    gene_tr_df = pd.read_table(gene2tr, header=None, names=["gene_id", "tr_id"])
    gene_interpro_df = gene_tr_df.merge(interpro_df)
    gene_go_df = gene_interpro_df[["gene_id", "go"]].dropna().copy()
    gene_go_df = gene_go_df.explode("go").drop_duplicates()
    gene_go_df.columns = ["gene", "term"]
    go_file = out_dir / "go.idmap.csv"
    gene_go_df.to_csv(go_file, index=False)

    gene_kegg_df = gene_interpro_df[["gene_id", "kegg"]].dropna().copy()
    gene_kegg_df = gene_kegg_df.explode("kegg").drop_duplicates()
    # gene_kegg_df.columns = ["ensembl", "ncbi"]
    gene_kegg_df.columns = ["gene", "term"]
    gene_kegg_df["ensembl"] = gene_kegg_df["gene"]
    gene_kegg_df["ncbi"] = gene_kegg_df["gene"]
    id_map_df = gene_kegg_df[["gene", "term"]].drop_duplicates()

    kegg_id_map = out_dir / "kegg.idmap.csv"
    kegg_gene_map = out_dir / "ncbi.idmap.csv"

    id_map_df.to_csv(kegg_id_map, index=False)
    ncbi_map_df = gene_kegg_df[["ensembl", "ncbi"]].drop_duplicates()
    ncbi_map_df.to_csv(kegg_gene_map, index=False, columns=["ensembl", "ncbi"])


if __name__ == "__main__":
    typer.run(main)
