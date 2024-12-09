import re
from pathlib import Path

import typer

MARKER_PATTERN = re.compile(r"\[(\[ATGC]+/[ATGC]+)\]")


from typing import Match


class SequenceError(Exception):
    """自定义序列错误异常类"""

    pass


class InvalidMarkerError(SequenceError):
    """无效标记错误"""

    pass


class InvalidBaseError(SequenceError):
    """无效碱基错误"""

    pass


# 定义SNP标记正则表达式模式
MARKER_PATTERN = re.compile(r"\[([ATGC]+/[ATGC]+)\]")


def validate_marker(sequence: str) -> Match:
    """
    验证序列中的SNP标记

    Args:
        sequence (str): 输入的DNA序列

    Returns:
        Match: 匹配的标记对象

    Raises:
        InvalidMarkerError: 当标记无效或未找到时
    """
    matches = list(MARKER_PATTERN.finditer(sequence))
    if not matches:
        raise InvalidMarkerError("未找到有效的SNP标记")
    if len(matches) > 1:
        raise SequenceError("序列中不允许包含多个SNP标记")
    return matches[0]


def validate_sequence_bases(sequence: str) -> None:
    """
    验证序列中的碱基

    Args:
        sequence (str): 输入的DNA序列

    Raises:
        InvalidBaseError: 当存在无效字符时
    """
    invalid_chars = set(sequence) - set("ATGC[]/")
    if invalid_chars:
        raise InvalidBaseError(f"序列包含无效字符: {', '.join(invalid_chars)}")


def extract_sequence(sequence: str) -> str:
    """
    从带有SNP标记的序列中提取参考序列

    Args:
        sequence (str): 输入的DNA序列，包含SNP标记 [REF/ALT]

    Returns:
        str: 提取的参考序列

    Raises:
        SequenceError: 当序列格式无效时

    Examples:
        >>> extract_sequence("ACGT[A/T]GCTA")
        'ACGTAGCTA'
    """
    # 验证序列基本格式
    if not sequence:
        raise SequenceError("序列不能为空")

    # 验证序列字符
    validate_sequence_bases(sequence)

    # 验证标记
    validate_marker(sequence)

    # 提取参考序列
    try:
        return MARKER_PATTERN.sub(lambda m: m.group(1).split("/")[0], sequence)
    except Exception as e:
        raise SequenceError(f"序列处理失败: {str(e)}") from e


def main(probe_fa: Path, out_fa: Path) -> None:
    with open(probe_fa, encoding="utf-8") as f:
        with open(out_fa, "w", encoding="utf-8") as out:
            for line in f:
                if line.startswith(">"):
                    out.write(f">{line.strip()}\n")
                else:
                    out.write(extract_sequence(line.strip()) + "\n")


if __name__ == "__main__":
    typer.run(main)
