from pathlib import Path

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


def main(out_dir: Path, chr_prefix: str = typer.Option(None)):
    for each_file in out_dir.glob("*assoc.txt"):
        df = pd.read_table(each_file)
        if chr_prefix:
            df["chr"] = df["chr"].map(lambda x: x.lstrip(chr_prefix))
        for each_test in TESTS:
            gt, model = each_file.stem.split(".")[:2]
            model_name = MODELS[model]
            test_name = TESTS[each_test]
            logger.info(f"Plotting {gt} {model_name} ({test_name})")
            manhattanplot_png = out_dir / f"{gt}_{model}_{each_test}.manhattan.png"
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
                manhattan_plot_prefix = out_dir / f"{gt}_{model}_{each_test}.manhattan"
                f.savefig(f"{manhattan_plot_prefix}.png", dpi=300)
                f.savefig(f"{manhattan_plot_prefix}.pdf")
                plt.close(f)
            qq_plot_png = out_dir / f"{gt}_{model}_{each_test}.QQ.png"
            if not qq_plot_png.is_file():
                ax = qqplot(
                    data=plot_df[f"p_{each_test}"],
                    title=f"{gt} {model_name} ({test_name})",
                )
                qq_plot_prefix = out_dir / f"{gt}_{model}_{each_test}.QQ"
                plt.savefig(f"{qq_plot_prefix}.png", dpi=300)
                plt.savefig(f"{qq_plot_prefix}.pdf")
                plt.close()


if __name__ == "__main__":
    typer.run(main)
