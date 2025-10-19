from pathlib import Path

import pandas as pd
import typer

CONTAINED_CLASS_CODE = ["=", "k", "c"]

app = typer.Typer()


def generate_class_code_df() -> pd.DataFrame:
    """生成 class_code 和 rank 对应表"""
    return pd.DataFrame(
        [{"class_code": code, "rank": n} for n, code in enumerate(CONTAINED_CLASS_CODE)]
    )


def load_gffcompare_tmap(file: Path) -> pd.DataFrame:
    """读取 gffcompare .tmap 文件的必要列"""
    return pd.read_table(
        file,
        usecols=[0, 2, 3],
        names=["qry_gene_id", "class_code", "ref_gene_id"],
        header=0,
    )


def process_overlap_mapping(
    gffcompare_df: pd.DataFrame, group_df: pd.DataFrame, name: str
) -> pd.DataFrame:
    """处理比对结果和 group 信息，输出映射表"""
    class_code_df = generate_class_code_df()

    overlap_df = class_code_df.merge(gffcompare_df, how="inner", on="class_code")
    overlap_df.sort_values(by=["rank"], inplace=True)
    overlap_df.drop_duplicates(subset="ref_gene_id", inplace=True)
    overlap_df.drop(columns=["rank", "class_code"], inplace=True)

    group_map_df = group_df[["group_id", "iwgsc_v2.1"]].rename(
        columns={"iwgsc_v2.1": "ref_gene_id"}
    )

    group_id_map_df = group_map_df.merge(overlap_df, on="ref_gene_id").drop(
        "ref_gene_id", axis=1
    )
    group_id_map_df.rename(columns={"qry_gene_id": name}, inplace=True)
    return group_id_map_df


@app.command()
def main(
    gffcompare_tmap_file: Path, group_info_file: Path, name: str, output_file: Path
):
    gffcompare_df = load_gffcompare_tmap(gffcompare_tmap_file)
    group_df = pd.read_table(group_info_file)
    result_df = process_overlap_mapping(gffcompare_df, group_df, name)
    result_df.to_csv(output_file, index=False)


if __name__ == "__main__":
    app()
    app()
