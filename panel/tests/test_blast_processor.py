from pathlib import Path

import pandas as pd
import pytest
import typer

from panel.blast_res_summary import (
    BlastConfig,
    calculate_match_statistics,
    filter_blast_results,
    get_blast_files,
    get_real_match_length,
    process_blast_file,
)


class TestBlastConfig:
    def test_default_values(self):
        """测试配置默认值"""
        config = BlastConfig()
        assert config.min_match_length == 24
        assert config.max_match_count == 1
        assert config.probe_length == 120
        assert not config.output_all
        assert not config.show_max_length
        assert config.min_identity is None
        assert config.max_mismatch_count is None
        assert config.max_gap_count is None

    def test_custom_values(self):
        """测试自定义配置值"""
        config = BlastConfig(
            min_match_length=30,
            max_match_count=2,
            probe_length=150,
            output_all=True,
            show_max_length=True,
            min_identity=95,
            max_mismatch_count=3,
            max_gap_count=1,
        )
        assert config.min_match_length == 30
        assert config.max_match_count == 2
        assert config.probe_length == 150
        assert config.output_all
        assert config.show_max_length
        assert config.min_identity == 95
        assert config.max_mismatch_count == 3
        assert config.max_gap_count == 1


class TestFileProcessing:
    def test_get_blast_files(self, temp_blast_dir):
        """测试获取 BLAST 文件列表"""
        files = get_blast_files(temp_blast_dir)
        assert len(files) == 1
        assert files[0].suffix == ".blast"

    def test_get_blast_files_empty_dir(self, tmp_path):
        """测试空目录情况"""
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        with pytest.raises(typer.Exit):
            get_blast_files(empty_dir)

    def test_get_blast_files_nonexistent_dir(self):
        """测试不存在的目录"""
        with pytest.raises(typer.Exit):
            get_blast_files(Path("/nonexistent/dir"))

    def test_process_blast_file(self, temp_blast_dir, blast_config):
        """测试处理单个 BLAST 文件"""
        blast_file = next(temp_blast_dir.glob("*.blast"))
        result = process_blast_file(blast_file, blast_config)
        assert isinstance(result, pd.DataFrame)
        assert not result.empty
        assert all(
            col in result.columns
            for col in ["id", "identity", "match_len", "mismatch", "gapopen"]
        )


class TestCalculations:
    def test_get_real_match_length(self):
        """测试实际匹配长度计算"""
        test_cases = [
            # match_len, mismatch, gapopen, probe_length, expected
            (100, 2, 0, 120, 98),  # 正常情况
            (150, 5, 1, 120, 116),  # 超过探针长度
            (50, 3, 1, 120, 46),  # 有错配和间隔
        ]

        for match_len, mismatch, gapopen, probe_length, expected in test_cases:
            row = pd.Series(
                {"match_len": match_len, "mismatch": mismatch, "gapopen": gapopen}
            )
            result = get_real_match_length(row, probe_length)
            assert result == expected

    def test_filter_blast_results(self, sample_blast_data, blast_config):
        """测试 BLAST 结果过滤"""
        filtered_df = filter_blast_results(sample_blast_data, blast_config)
        assert len(filtered_df) <= len(sample_blast_data)
        assert all(filtered_df["match_len"] >= blast_config.min_match_length)
        if blast_config.min_identity:
            assert all(filtered_df["identity"] >= blast_config.min_identity)

    def test_calculate_match_statistics(self, sample_blast_data, blast_config):
        """测试匹配统计计算"""
        result = calculate_match_statistics(sample_blast_data, blast_config)

        assert isinstance(result, pd.DataFrame)
        assert "id" in result.columns
        assert "blast_match" in result.columns
        assert "second_max_length" in result.columns

        # 测试空数据框情况
        empty_df = pd.DataFrame(columns=sample_blast_data.columns)
        empty_result = calculate_match_statistics(empty_df, blast_config)
        assert len(empty_result) == 0
        assert all(
            col in empty_result.columns
            for col in ["id", "blast_match", "second_max_length"]
        )

    def test_calculate_match_statistics_with_max_length(
        self, sample_blast_data, blast_config
    ):
        """测试带最大长度的匹配统计计算"""
        blast_config.show_max_length = True
        result = calculate_match_statistics(sample_blast_data, blast_config)

        assert "max_length" in result.columns
        assert all(result["max_length"] >= result["second_max_length"])
