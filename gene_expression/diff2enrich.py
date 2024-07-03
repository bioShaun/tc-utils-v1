from pathlib import Path

import delegator
import pandas as pd
import typer
from loguru import logger

ENRICH_PY = Path(__file__).parent / "enrichment.py"


def filter_diff_genes(
    diff_df: pd.DataFrame,
    compare: str,
    out_dir: Path,
    logFc: float = 1,
    pval: float = 0.05,
    is_adjust_p: bool = True,
) -> None:
    pval_col = "FDR" if is_adjust_p else "PValue"
    mask1 = diff_df[pval_col] <= pval
    mask2 = diff_df["logFC"] >= logFc
    mask3 = diff_df["logFC"] <= -logFc
    up_diff = diff_df[mask1 & mask2]
    down_diff = diff_df[mask1 & mask3]
    all_diff = pd.concat([up_diff, down_diff])
    diff_dir = out_dir / "diff"
    diff_dir.mkdir(parents=True, exist_ok=True)
    up_diff_file = diff_dir / f"{compare}.UP.diff_genes.txt"
    down_diff_file = diff_dir / f"{compare}.DOWN.diff_genes.txt"
    all_diff_file = diff_dir / f"{compare}.ALL.diff_genes.txt"
    up_diff.to_csv(
        up_diff_file, sep="\t", index=False, header=None, columns=["Gene_ID"]
    )
    down_diff.to_csv(
        down_diff_file, sep="\t", index=False, header=None, columns=["Gene_ID"]
    )
    all_diff.to_csv(
        all_diff_file, sep="\t", index=False, header=None, columns=["Gene_ID"]
    )


def go_enrichment(
    diff_gene_list: Path, out_prefix: Path, go_id_map: Path, go_des: Path
) -> None:
    cmd = f"python {ENRICH_PY} go {diff_gene_list}  {go_id_map} {go_des} {out_prefix}"
    print(cmd)
    delegator.run(cmd)


def kegg_enrichment(
    diff_gene_list: Path,
    out_prefix: Path,
    kegg_id_map: Path,
    kegg_des: Path,
    kegg_abbr: str,
    ncbi_id_map: Path,
):
    cmd = f"python {ENRICH_PY} kegg {diff_gene_list} {kegg_id_map} {kegg_des} {out_prefix} {kegg_abbr} --ncbi-map {ncbi_id_map}"
    print(cmd)
    delegator.run(cmd)


def main(
    diff_dir: Path,
    out_dir: Path,
    go_id_map: Path,
    go_des: Path,
    kegg_id_map: Path,
    kegg_des: Path,
    kegg_abbr: str,
    ncbi_id_map: Path,
    logFc: float = 1,
    pval: float = 0.05,
    is_adjust_p: bool = True,
    run_go: bool = True,
    run_kegg: bool = True,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for diff_dir_i in diff_dir.iterdir():
        if not diff_dir_i.is_dir():
            continue
        compare_name = diff_dir_i.name
        diff_file = diff_dir_i / f"{compare_name}.edgeR.DE_results.txt"
        if not diff_file.exists():
            raise ValueError(f"{diff_file} not exists")
        diff_df = pd.read_csv(diff_file, sep="\t")
        logger.info(f"run {compare_name}")
        filter_diff_genes(diff_df, compare_name, out_dir, logFc, pval, is_adjust_p)
        for reg in ["ALL", "UP", "DOWN"]:
            diff_dir = out_dir / "diff"
            go_enrich_dir = out_dir / "enrichment" / "go" / compare_name
            kegg_enrich_dir = out_dir / "enrichment" / "kegg" / compare_name
            go_enrich_dir.mkdir(parents=True, exist_ok=True)
            kegg_enrich_dir.mkdir(parents=True, exist_ok=True)
            diff_gene_list = diff_dir / f"{compare_name}.{reg}.diff_genes.txt"
            go_out = go_enrich_dir / f"{compare_name}.{reg}.go_enrichment"
            kegg_out = kegg_enrich_dir / f"{compare_name}.{reg}.kegg_enrichment"
            if run_go:
                go_enrichment(
                    diff_gene_list,
                    go_out,
                    go_id_map,
                    go_des,
                )
            if run_kegg:
                kegg_enrichment(
                    diff_gene_list,
                    kegg_out,
                    kegg_id_map,
                    kegg_des,
                    kegg_abbr,
                    ncbi_id_map,
                )


if __name__ == "__main__":
    typer.run(main)
