from pathlib import Path

import pandas as pd
import typer

OVERLAP_CLASS_CODE_ALL = ["=", "k", "c", "m", "j", "o"]
CONTAINED_CLASS_CODE = ["=", "k", "c"]


def generate_class_code_df() -> pd.DataFrame:
    class_code_dict_list = []
    for n, code_i in enumerate(CONTAINED_CLASS_CODE):
        class_code_dict_list.append({"class_code": code_i, "rank": n})
    return pd.DataFrame(class_code_dict_list)


def main(
    gffcompare_tmap_file: Path, group_info_file: Path, name: str, output_file: Path
):
    class_code_df = generate_class_code_df()
    gff_compare_df = pd.read_table(gffcompare_tmap_file, usecols=[0, 2, 3])
    overlap_df = class_code_df.merge(gff_compare_df, how="inner", on="class_code")
    overlap_df.sort_values(by=["rank"], inplace=True)
    overlap_df.drop_duplicates(subset="ref_gene_id", inplace=True)
    overlap_df.drop(columns=["rank", "class_code"], inplace=True)
    group_df = pd.read_table(group_info_file)
    group_map_df = group_df[["group_id", "iwgsc_v2.1"]].copy()
    group_map_df.rename(columns={"iwgsc_v2.1": "ref_gene_id"}, inplace=True)
    group_id_map_df = group_map_df.merge(overlap_df, on="ref_gene_id").drop(
        "ref_gene_id", axis=1
    )
    group_id_map_df.rename(columns={"qry_gene_id": name}, inplace=True)
    group_id_map_df.to_csv(output_file, index=False)


if __name__ == "__main__":
    typer.run(main)
