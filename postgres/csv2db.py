import pandas as pd
import psycopg2
import typer

from pathlib import Path
from sqlalchemy import create_engine


def formatTitle(title) -> str:
    return title.replace("-", "_").lower()


def toPsql(
    csv: Path,
    name: str,
    username: str,
    password: str,
    db: str,
    host="localhost",
    port=5432,
    add_id: bool = typer.Option(False, "--add-id"),
    replace_table: bool = typer.Option(False, "--replace"),
    format_title: bool = typer.Option(False, "--format-title"),
    sep=",",
):
    engine = create_engine(
        f"postgresql+psycopg2://{username}:{password}:@{host}:{port}/{db}"
    )
    df = pd.read_csv(csv, sep=sep)
    if format_title:
        df.columns = [formatTitle(each) for each in df.columns]
    if add_id:
        df.loc[:, "id"] = [each + 1 for each in df.index]
    if replace_table:
        df.to_sql(name, engine, if_exists="replace", index=False, chunksize=10000)
    else:
        df.to_sql(name, engine, if_exists="append", index=False, chunksize=10000)


if __name__ == "__main__":
    typer.run(toPsql)
