from pathlib import Path

import pandas as pd
import typer


def main(gt: Path, id_map: Path, out: Path):
    gt_df = pd.read_table(gt)
    id_map_df = pd.read_table(id_map, header=None, names=["probe_id", "id"])
    gt_df["id"] = gt_df["CHROM"].str.cat(gt_df["POS"].astype(str), sep="_")
    probe_gt_df = id_map_df.merge(gt_df).drop("id", axis=1)
    probe_gt_df.to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)
