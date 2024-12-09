import re
from pathlib import Path
from typing import Match, Set

import typer


class SequenceError(Exception):
    """自定义序列错误异常类"""

    pass


class InvalidMarkerError(SequenceError):
    """无效标记错误"""

    pass


class InvalidBaseError(SequenceError):
    """无效碱基错误"""

    pass


# IUPAC核苷酸编码对照表
IUPAC_CODES = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "M": {"A", "C"},
    "K": {"G", "T"},
    "S": {"C", "G"},
    "W": {"A", "T"},
    "H": {"A", "C", "T"},
    "B": {"C", "G", "T"},
    "V": {"A", "C", "G"},
    "D": {"A", "G", "T"},
    "N": {"A", "C", "G", "T"},
}

# 定义SNP标记正则表达式模式
MARKER_PATTERN = re.compile(r"\[([A-Za-z]+/[A-Za-z]+)\]")


def get_iupac_bases(iupac_code: str) -> Set[str]:
    """获取IUPAC码对应的碱基集合"""
    try:
        return IUPAC_CODES[iupac_code.upper()]
    except KeyError:
        raise InvalidBaseError(f"无效的IUPAC码: {iupac_code}")


def validate_marker(sequence: str) -> Match:
    """
    验证序列中的SNP标记
    
    Args:
        sequence (str): 输入的DNA序列
        
    Returns:
        Match: 匹配的标记对象
        
    Raises:
        InvalidMarkerError: 当标记无效或未找到时
        SequenceError: 当存在多个标记时
    """
    matches = list(MARKER_PATTERN.finditer(sequence))
    if not matches:
        raise InvalidMarkerError("未找到有效的SNP标记")
    if len(matches) > 1:
        raise SequenceError("序列中不允许包含多个SNP标记")
    
    # 验证SNP标记中的碱基
    match = matches[0]
    alleles = match.group(1).split('/')
    if len(alleles) != 2:
        raise InvalidMarkerError("SNP标记格式错误")
    
    # 验证每个碱基是否为有效的IUPAC码
    for allele in alleles:
        for base in allele:
            if base.upper() not in IUPAC_CODES:
                raise InvalidBaseError(f"无效的IUPAC码: {base}")
    
    return match


def validate_sequence_bases(sequence: str, marker_match: Match) -> None:
    """验证序列中的碱基"""
    # 移除SNP标记部分
    start, end = marker_match.span()
    sequence_without_marker = sequence[:start] + sequence[end:]

    # 检查剩余序列中的碱基
    invalid_chars = set(sequence_without_marker) - set(IUPAC_CODES.keys())
    if invalid_chars:
        raise InvalidBaseError(f"序列包含无效字符: {', '.join(sorted(invalid_chars))}")


def extract_sequence(sequence: str) -> str:
    """从带有SNP标记的序列中提取参考序列"""
    if not sequence:
        raise SequenceError("序列不能为空")

    # 先验证标记
    match = validate_marker(sequence)
    # 再验证其他碱基
    validate_sequence_bases(sequence, match)

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
