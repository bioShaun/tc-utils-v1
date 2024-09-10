from pathlib import Path

import pandas as pd
import psycopg2
import typer
from sqlalchemy import create_engine
from tqdm import tqdm


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
    format_title: bool = typer.Option(False, "--format-title"),
    sep=",",
    chunk_size: int = 10000,
):
    engine = create_engine(
        f"postgresql+psycopg2://{username}:{password}:@{host}:{port}/{db}"
    )
    dfs = pd.read_csv(csv, sep=sep, chunksize=chunk_size)
    for df in tqdm(dfs):
        if format_title:
            df.columns = [formatTitle(each) for each in df.columns]
        if add_id:
            df.loc[:, "id"] = [each + 1 for each in df.index]
        else:
            df.to_sql(name, engine, if_exists="append", index=False, chunksize=10000)


if __name__ == "__main__":
    typer.run(toPsql)
