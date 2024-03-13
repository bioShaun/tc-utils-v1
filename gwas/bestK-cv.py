from pathlib import Path

import pandas as pd
import seaborn as sns
import typer


def k_table(structure_dir: Path) -> pd.DataFrame:
    """
    Create a table of likelihoods for different values of K.
    """
    dict_list = []
    for structure_path in structure_dir.glob("*.log"):
        k_value = int(structure_path.name.split(".")[2])
        with structure_path.open("r") as f:
            for eachline in f:
                if "CV error" in eachline:
                    cv = float(eachline.split(' = ')[-1].split(',')[0])
                    #cv = float(eachline.split()[-1])
                    dict_list.append({"K": k_value, "cv": cv})
    return pd.DataFrame(dict_list)


def main(structure_dir: Path, out_path: Path) -> None:
    """
    Plot the likelihood of different values of K.
    """
    df = k_table(structure_dir=structure_dir)
    bestK = df.loc[df["cv"].idxmin(), "K"]
    fig_length = len(df) * 0.75
    g = sns.lineplot(data=df, x="K", y="cv", marker="o")

    g.set(xticks=range(min(df["K"]), max(df["K"]) + 1, 1))  # type: ignore
    g.set_ylabel("CV error")
    ymin, ymax = g.get_ylim()
    g.vlines(
        x=[bestK],
        ymin=ymin,
        ymax=ymax,
        colors=["tab:orange", "tab:blue"],
        ls="--",
        lw=2,
    )
    fig = g.get_figure()
    fig.set_figwidth(fig_length)
    fig.set_figheight(6)
    out_path.mkdir(exist_ok=True, parents=True)
    fig.savefig(f"{out_path}/cv-error.png")


if __name__ == "__main__":
    typer.run(main)
