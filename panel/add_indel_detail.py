import re
from enum import Enum
from pathlib import Path

import pandas as pd
import pytest
import typer


class IndelType(Enum):
    SNP = "SNP"
    DEL = "DEL"
    INS = "INS"


def get_ref_alt_lengths(alleles: str) -> tuple[int, list[int]]:
    """
    Given a string of alleles, return the length of the reference allele and a list of
    lengths of the alternative alleles.

    Args:
        alleles: A string of alleles, separated by a forward slash and commas.

    Returns:
        A tuple where the first element is the length of the reference allele and the
        second element is a list of lengths of the alternative alleles.
    """
    alleles = re.sub("-|del", "", alleles)
    if not alleles or "/" not in alleles:
        raise ValueError("Invalid alleles format")
    ref, alt = alleles.split("/")

    if not ref and not alt:
        raise ValueError("Both reference and alternative alleles are empty")
    ref_len = len(ref)
    alt_lens = [len(a) for a in alt.split(",")]
    return ref_len, alt_lens


def get_indel_type(alleles: str) -> IndelType:
    """
    Determines the type of indel (insertion, deletion, or SNP) based on the provided alleles string.

    Args:
        alleles (str): A string representing the reference and alternative alleles, separated by a forward slash (e.g. "A/T,C").

    Returns:
        IndelType: An enum representing the type of indel, either IndelType.SNP, IndelType.DEL, or IndelType.INS.
    """

    ref_len, alt_lens = get_ref_alt_lengths(alleles)

    # Determine if all alternative alleles are the same length as the reference
    if all(ref_len == alt_len for alt_len in alt_lens):
        return IndelType.SNP

    # Check if any alternative allele is shorter than the reference
    if any(ref_len > alt_len for alt_len in alt_lens):
        return IndelType.DEL

    # Otherwise, it's an insertion
    return IndelType.INS


def get_indel_length(alleles: str) -> int:
    """
    Calculates the length of an indel given the alleles string.

    The length of an indel is the maximum difference in length between the reference and
    alternative alleles.

    Args:
        alleles (str): A string representing the reference and alternative alleles,
            separated by a forward slash (e.g. "A/T,C").

    Returns:
        int: The length of the indel.
    """
    ref_len, alt_lengths = get_ref_alt_lengths(alleles)
    max_delta = max(abs(ref_len - alt_len) for alt_len in alt_lengths)
    return max_delta


def get_row_indel_length(row: pd.Series) -> int:
    return get_indel_length(row["alleles"])


@pytest.mark.parametrize(
    "alleles, expected",
    [
        ("ATG/ATGC,ATGCA", 2),
        ("ATG/AT", 1),
        ("ATG/A,TG,ATGCCC", 3),
        ("ATG/ATG", 0),  # 没有插入或缺失
    ],
)
def test_get_indel_length(alleles, expected):
    assert get_indel_length(alleles) == expected


# 测试代码
@pytest.mark.parametrize(
    "alleles, expected_type",
    [
        ("A/T", IndelType.SNP),
        ("G/C", IndelType.SNP),
        ("AT/AT", IndelType.SNP),
        ("A/A,T", IndelType.SNP),
    ],
)
def test_snp(alleles, expected_type):
    assert get_indel_type(alleles) == expected_type


@pytest.mark.parametrize(
    "alleles, expected_type",
    [
        ("A/-", IndelType.DEL),
        ("A/del", IndelType.DEL),
        ("AT/A", IndelType.DEL),
        ("ATG/A", IndelType.DEL),
        ("ATGC/AT", IndelType.DEL),
        ("AT/A,ATG", IndelType.DEL),
    ],
)
def test_deletion(alleles, expected_type):
    assert get_indel_type(alleles) == expected_type


@pytest.mark.parametrize(
    "alleles, expected_type",
    [
        ("A/AT", IndelType.INS),
        ("G/GTT", IndelType.INS),
        ("AT/ATG", IndelType.INS),
        ("A/AT,ATG", IndelType.INS),
    ],
)
def test_insertion(alleles, expected_type):
    assert get_indel_type(alleles) == expected_type


def test_complex_cases():
    assert get_indel_type("A/T,AT,ATG") == IndelType.INS
    assert get_indel_type("ATG/A,AT,ATGC") == IndelType.DEL
    assert get_indel_type("AT/A,ATG,AT") == IndelType.DEL


def test_edge_cases():
    with pytest.raises(ValueError):
        get_ref_alt_lengths("")
    with pytest.raises(ValueError):
        get_ref_alt_lengths("A")
    with pytest.raises(ValueError):
        get_ref_alt_lengths("/")


def get_row_indel_type(row: pd.Series) -> IndelType:
    return get_indel_type(row["alleles"])


def main(pos_file: Path, out_file: Path):
    df = pd.read_table(pos_file)
    df["indel_type"] = df.apply(get_row_indel_type, axis=1)  # type: ignore
    df["indel_length"] = df.apply(get_row_indel_length, axis=1)  # type: ignore
    df.to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)
