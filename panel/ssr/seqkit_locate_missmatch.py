import pandas as pd
import typer


def get_mismatch_count(seq: str, pattern: str) -> int:
    """
    计算字符串中模式字符的个数

    参数:
    seq -- 待匹配的字符串
    pattern -- 模式字符

    返回:
    模式字符的个数
    """
    if len(pattern) != len(seq):
        raise ValueError("序列长度不同，无法逐位比较。")
    return sum(a != b for a, b in zip(pattern.upper(), seq.upper()))


def main(
    seqkit_locate_table: Path = typer.Argument(..., help="seqkit locate table")
) -> None:
    df = pd.read_table(seqkit_locate_table)
    df["mismatch"] = df.apply(
        lambda x: get_mismatch_count(x["pattern"], x["matched"]), axis=1
    )
    out_file = seqkit_locate_table.with_suffix(".add-mismatch.tsv")
    df.to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)
