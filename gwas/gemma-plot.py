from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd
import typer
from loguru import logger
from qmplot import manhattanplot, qqplot

TESTS = {
    "wald": "Wald test",
    "lrt": "Likelihood ratio test",
    "score": "Score test",
}

MODELS = {
    "lm": "Linear model",
    "lmm": "Linear mixed model",
}


def main(
    data_dir: Path,
    out_dir: Optional[Path] = None,
    chr_prefix: str = typer.Option(None),
    chr_list_file: Optional[Path] = None,
):
    for each_file in data_dir.glob("*/*assoc.txt"):
        file_path = each_file.parent
        df = pd.read_table(each_file)
        df["chr"] = df["chr"].astype(str)
        if chr_list_file:
            chr_list = pd.read_table(chr_list_file, header=None)[0].to_list()
            df = df[df["chr"].isin(chr_list)]
        plot_dir = file_path if out_dir is None else out_dir
        if chr_prefix:
            df["chr"] = df["chr"].map(lambda x: x.lstrip(chr_prefix))
        for each_test in TESTS:
            gt, model = each_file.stem.split(".")[:2]
            model_name = MODELS[model]
            test_name = TESTS[each_test]
            logger.info(f"Plotting {gt} {model_name} ({test_name})")
            manhattanplot_png = plot_dir / f"{gt}_{model}_{each_test}.manhattan.png"
            plot_df = df.dropna(subset=[f"p_{each_test}"]).copy()
            if not manhattanplot_png.is_file():
                f, ax = plt.subplots(figsize=(12, 4), facecolor="w", edgecolor="k")
                manhattanplot(
                    data=plot_df,
                    chrom="chr",
                    pos="ps",
                    pv=f"p_{each_test}",
                    ax=ax,
                    title=f"{gt} {model_name} ({test_name})",
                    hline_kws={"linestyle": "--", "lw": 1.3},
                    sign_line_cols=["#D62728", "#2CA02C"],
                )
                ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
                manhattan_plot_prefix = plot_dir / f"{gt}_{model}_{each_test}.manhattan"
                plt.tight_layout()
                # f.subplots_adjust(bottom=0.25)
                f.savefig(f"{manhattan_plot_prefix}.png", dpi=300)
                f.savefig(f"{manhattan_plot_prefix}.pdf")
                plt.close(f)
            qq_plot_png = plot_dir / f"{gt}_{model}_{each_test}.QQ.png"
            if not qq_plot_png.is_file():
                plt.figure(figsize=(12, 6))
                ax = qqplot(
                    data=plot_df[f"p_{each_test}"],
                    title=f"{gt} {model_name} ({test_name})",
                )
                ax.title.set_fontsize(8)
                qq_plot_prefix = plot_dir / f"{gt}_{model}_{each_test}.QQ"
                plt.savefig(f"{qq_plot_prefix}.png", dpi=300)
                plt.savefig(f"{qq_plot_prefix}.pdf")
                plt.close()


if __name__ == "__main__":
    typer.run(main)
