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
):
    engine = create_engine(
        f"postgresql+psycopg2://{username}:{password}:@{host}:{port}/{db}"
    )
    df = pd.read_sql(name, engine)
    df.to_csv(csv, index=False)


if __name__ == "__main__":
    typer.run(toPsql)
