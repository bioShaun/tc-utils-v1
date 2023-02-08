import pandas as pd
import psycopg2
import typer

from pathlib import Path
from sqlalchemy import create_engine


def toPsql(
    csv: Path,
    name: str,
    username: str,
    password: str,
    db: str,
    host="localhost",
    port=5432,
    add_id=False,
    append=False,
    sep=",",
):
    engine = create_engine(
        f"postgresql+psycopg2://{username}:{password}:@{host}:{port}/{db}"
    )
    df = pd.read_csv(csv, sep=sep)
    if add_id:
        df.loc[:, "id"] = [each + 1 for each in df.index]
    if append:
        df.to_sql(name, engine, if_exists="append", index=False, chunksize=10000)
    else:
        df.to_sql(name, engine, if_exists="replace", index=False, chunksize=10000)


if __name__ == "__main__":
    typer.run(toPsql)
