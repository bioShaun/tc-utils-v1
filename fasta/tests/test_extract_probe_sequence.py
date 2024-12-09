import pytest

from fasta.extract_probe_sequence import (
    InvalidBaseError,
    InvalidMarkerError,
    SequenceError,
    extract_sequence,
)


class TestSequenceExtractor:
    """序列提取器测试类"""

    def test_valid_sequence(self):
        """测试有效序列"""
        test_cases = [
            {"input": "ACGT[A/T]GCTA", "expected": "ACGTAGCTA"},
            {"input": "GAGTTC[C/A]ATGGAG", "expected": "GAGTTCCATGGAG"},
            {"input": "[A/G]TCGA", "expected": "ATCGA"},
            {"input": "TCGA[G/C]", "expected": "TCGAG"},
        ]

        for case in test_cases:
            assert extract_sequence(case["input"]) == case["expected"]

    def test_empty_sequence(self):
        """测试空序列"""
        with pytest.raises(SequenceError) as exc_info:
            extract_sequence("")
        assert str(exc_info.value) == "序列不能为空"

    def test_invalid_bases(self):
        """测试无效碱基"""
        test_cases = [
            "ACGT$[A/T]GCTA",
            "ACGT[A/T]GCT@",
            "AC1GT[A/T]GCTA",
        ]

        for sequence in test_cases:
            with pytest.raises(InvalidBaseError):
                extract_sequence(sequence)

    def test_invalid_marker(self):
        """测试无效标记"""
        test_cases = [
            "ACGTGCTA",  # 无标记
            "ACGT[X/Y]GCTA",  # 非法碱基标记
            "ACGT[AT]GCTA",  # 错误的标记格式
        ]

        for sequence in test_cases:
            with pytest.raises((InvalidMarkerError, InvalidBaseError)):
                extract_sequence(sequence)

    def test_multiple_markers(self):
        """测试多个标记"""
        sequence = "ACGT[A/T]GCTA[C/G]TGCA"
        with pytest.raises(SequenceError):
            extract_sequence(sequence)

    @pytest.mark.parametrize(
        "sequence,expected",
        [
            ("ACGT[A/T]GCTA", "ACGTAGCTA"),
            ("GAGTTC[C/A]ATGGAG", "GAGTTCCATGGAG"),
            ("[A/G]TCGA", "ATCGA"),
        ],
    )
    def test_parametrize_valid_sequences(self, sequence, expected):
        """参数化测试有效序列"""
        assert extract_sequence(sequence) == expected
