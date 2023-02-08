from pathlib import Path
import delegator
import typer
import pandas as pd

from io import StringIO

app = typer.Typer()


BEDTOOLS_OUT_COLUMNS = [
    "w_chrom",
    "w_start",
    "w_end",
    "count",
    "chrom",
    "start",
    "end",
    "ref",
    "alt",
    "num1",
    "num2",
    "overlap",
]

OUT_COLUMNS = ["chrom", "start", "end", "ref", "alt", "num1", "num2"]


@app.command()
def random_select(var_table: Path, region_table: Path, out: Path, seed: int = 1):
    bedtools_cmd = f'awk \'{{print $1"\\t"$2-1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6}}\' {var_table} |bedtools intersect -a {region_table} -b - -wo'
    res = delegator.run(bedtools_cmd)
    if res.err:
        print(res.err)
    if res.out:    
        df = pd.read_csv(
            StringIO(res.out), sep="\t", header=None, names=BEDTOOLS_OUT_COLUMNS
        )
        random_df = df.groupby(["w_chrom", "w_start"]).sample(n=1, random_state=seed)
        random_df.to_csv(out, sep="\t", header=False, columns=OUT_COLUMNS, index=False)


@app.command()
def select():
    print("select")


if __name__ == "__main__":
    app()
