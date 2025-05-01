import pytest

from fasta.extract_probe_sequence import (
    InvalidBaseError,
    SequenceError,
    extract_sequence,
    get_iupac_bases,
)


class TestSequenceExtractor:
    """序列提取器测试类"""

    def test_valid_sequence_with_iupac(self):
        """测试包含IUPAC码的有效序列"""
        test_cases = [
            {"input": "ACGT[R/Y]GCTA", "expected": "ACGTRGCTA"},
            {"input": "GAGTTC[N/A]ATGGAG", "expected": "GAGTTCNATGGAG"},
            {"input": "[M/K]TCGA", "expected": "MTCGA"},
            {"input": "TCGA[S/W]", "expected": "TCGAS"},
        ]

        for case in test_cases:
            result = extract_sequence(case["input"])
            assert (
                result == case["expected"]
            ), f"输入: {case['input']}, 期望: {case['expected']}, 实际: {result}"

    def test_invalid_iupac_code(self):
        """测试无效的IUPAC码"""
        test_cases = [
            "ACGT[X/T]GCTA",
            "ACGT[A/Z]GCTA",
            "ACGT[P/Q]GCTA",
        ]

        for sequence in test_cases:
            with pytest.raises(InvalidBaseError, match="无效的IUPAC码"):
                extract_sequence(sequence)

    def test_get_iupac_bases(self):
        """测试IUPAC码转换"""
        test_cases = [
            ("A", {"A"}),
            ("R", {"A", "G"}),
            ("N", {"A", "C", "G", "T"}),
            ("Y", {"C", "T"}),
        ]

        for iupac_code, expected_bases in test_cases:
            assert get_iupac_bases(iupac_code) == expected_bases

    def test_invalid_iupac_lookup(self):
        """测试查找无效的IUPAC码"""
        with pytest.raises(InvalidBaseError, match="无效的IUPAC码"):
            get_iupac_bases("X")

    # 保留之前的其他测试方法...
    def test_empty_sequence(self):
        """测试空序列"""
        with pytest.raises(SequenceError, match="序列不能为空"):
            extract_sequence("")

    def test_multiple_markers(self):
        """测试多个标记"""
        sequence = "ACGT[R/Y]GCTA[M/K]TGCA"
        with pytest.raises(SequenceError, match="序列中不允许包含多个SNP标记"):
            extract_sequence(sequence)

    @pytest.mark.parametrize(
        "sequence,expected",
        [
            ("ACGT[R/Y]GCTA", "ACGTRGCTA"),
            ("GAGTTC[N/A]ATGGAG", "GAGTTCNATGGAG"),
            ("[M/K]TCGA", "MTCGA"),
        ],
    )
    def test_parametrize_valid_sequences(self, sequence, expected):
        """参数化测试有效序列"""
        assert extract_sequence(sequence) == expected
