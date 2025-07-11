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

    if alleles.count("/") > 1:
        raise ValueError("Invalid alleles format, more than one '/'")

    ref_len = len(ref)
    alt_lens = [len(a) for a in alt.split(",")]
    return ref_len, alt_lens


def get_indel_type(alleles: str) -> str:
    """
    Given a string of alleles, determine the type of indel.

    Args:
        alleles: A string of alleles, separated by a forward slash and commas.

    Returns:
        A string indicating the type of indel, one of "SNP", "DEL" or "INS".
    """
    ref_len, alt_lens = get_ref_alt_lengths(alleles)

    # Determine if all alternative alleles are the same length as the reference
    if all(ref_len == alt_len for alt_len in alt_lens):
        return IndelType.SNP.value

    # Check if any alternative allele is shorter than the reference
    if any(ref_len > alt_len for alt_len in alt_lens):
        return IndelType.DEL.value

    # Otherwise, it's an insertion
    return IndelType.INS.value


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
        ("A/T", IndelType.SNP.value),
        ("G/C", IndelType.SNP.value),
        ("AT/AT", IndelType.SNP.value),
        ("A/A,T", IndelType.SNP.value),
    ],
)
def test_snp(alleles, expected_type):
    assert get_indel_type(alleles) == expected_type


@pytest.mark.parametrize(
    "alleles, expected_type",
    [
        ("A/-", IndelType.DEL.value),
        ("A/del", IndelType.DEL.value),
        ("AT/A", IndelType.DEL.value),
        ("ATG/A", IndelType.DEL.value),
        ("ATGC/AT", IndelType.DEL.value),
        ("AT/A,ATG", IndelType.DEL.value),
    ],
)
def test_deletion(alleles, expected_type):
    assert get_indel_type(alleles) == expected_type


@pytest.mark.parametrize(
    "alleles, expected_type",
    [
        ("A/AT", IndelType.INS.value),
        ("G/GTT", IndelType.INS.value),
        ("AT/ATG", IndelType.INS.value),
        ("A/AT,ATG", IndelType.INS.value),
    ],
)
def test_insertion(alleles, expected_type):
    assert get_indel_type(alleles) == expected_type


def test_complex_cases():
    assert get_indel_type("A/T,AT,ATG") == IndelType.INS.value
    assert get_indel_type("ATG/A,AT,ATGC") == IndelType.DEL.value
    assert get_indel_type("AT/A,ATG,AT") == IndelType.DEL.value


def test_edge_cases():
    with pytest.raises(ValueError):
        get_ref_alt_lengths("")
    with pytest.raises(ValueError):
        get_ref_alt_lengths("A")
    with pytest.raises(ValueError):
        get_ref_alt_lengths("/")
    with pytest.raises(ValueError):
        get_ref_alt_lengths("A/T/GG")


def get_row_indel_type(row: pd.Series) -> str:
    return get_indel_type(row["alleles"])


def main(pos_file: Path):
    df = pd.read_table(pos_file)
    df["ref"] = df["alleles"].map(lambda x: x.split("/")[0])
    df["ref_length"] = df["ref"].map(len)
    center_df = df[df["ref_length"] == 1]
    center_file = pos_file.with_suffix(".center.tsv")
    center_df.to_csv(center_file, sep="\t", index=False)
    side_df = df[df["ref_length"] > 1]
    right_df = side_df.copy()
    right_df["id"] = right_df.apply(lambda x: f"{x['chrom']}_{x['pos']}", axis=1)
    left_file = pos_file.with_suffix(".left.tsv")
    left_df = side_df.copy()
    left_df.drop(["ref", "ref_length"], axis=1, inplace=True)
    left_df.to_csv(left_file, sep="\t", index=False)

    right_df["pos"] = right_df["pos"] + right_df["ref_length"] - 1
    right_df.drop(["ref", "ref_length"], axis=1, inplace=True)
    right_file = pos_file.with_suffix(".right.tsv")
    right_df.to_csv(right_file, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)
