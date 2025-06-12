#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, Mock, mock_open, patch

import numpy as np
import pandas as pd
import pytest

# 假设主模块名为 fastq_processor
from fq.project_fq import (
    FASTQ_EXTENSIONS,
    READ_TYPE_PATTERNS,
    FastqProcessor,
    ScriptRunner,
    log_statistics,
    validate_sample_info,
    write_nextflow_input,
)


class TestFastqProcessor:
    """测试FastqProcessor类"""

    def setup_method(self):
        """每个测试方法前的设置"""
        self.temp_dir = Path(tempfile.mkdtemp())
        self.processor = FastqProcessor(self.temp_dir)

    def teardown_method(self):
        """每个测试方法后的清理"""
        import shutil

        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_init_valid_directory(self):
        """测试有效目录初始化"""
        assert self.processor.base_dir == self.temp_dir
        assert self.processor.base_dir.exists()

    def test_init_invalid_directory(self):
        """测试无效目录初始化"""
        invalid_dir = Path("/nonexistent/directory")
        with pytest.raises(FileNotFoundError, match="基础目录不存在"):
            FastqProcessor(invalid_dir)

    def test_determine_read_type(self):
        """测试读取类型识别"""
        test_cases = [
            ("sample_combined_R1.fastq.gz", "R1"),
            ("sample_combined_R2.fastq.gz", "R2"),
            ("sample_R1.fastq.gz", "R1"),
            ("sample_R2.fastq.gz", "R2"),
            ("sample_1.fastq.gz", "R1"),
            ("sample_2.fastq.gz", "R2"),
            ("sample_unknown.fastq.gz", None),
        ]

        for filename, expected in test_cases:
            result = self.processor._determine_read_type(filename)
            assert result == expected, f"Failed for {filename}"

    @patch("pathlib.Path.exists", return_value=True)
    @patch("pathlib.Path.glob")
    def test_parse_fastq_filename_fq_extension(self, mock_glob, mock_exists):
        processor = FastqProcessor(Path("/fake/base"))

        # 创建 FASTQ 文件 mock
        mock_fastq1 = Mock()
        mock_fastq1.name = "sample_1.fq.gz"
        mock_fastq1.absolute.return_value = Path("/path/to/sample_1.fq.gz")

        mock_fastq2 = Mock()
        mock_fastq2.name = "sample_2.fq.gz"
        mock_fastq2.absolute.return_value = Path("/path/to/sample_2.fq.gz")

        # 设置 glob 的 side effect
        def glob_side_effect(pattern):
            if pattern == "*.fastq.gz":
                return []
            elif pattern == "*.fq.gz":
                return [mock_fastq1, mock_fastq2]
            return []

        mock_glob.side_effect = glob_side_effect

        sample_path = Path("Sample-LIB002")
        result = processor.parse_fastq_filename(sample_path)

        expected = [
            {
                "libid": "LIB002",
                "read_type": "R1",
                "path": "/path/to/sample_1.fq.gz",
            },
            {
                "libid": "LIB002",
                "read_type": "R2",
                "path": "/path/to/sample_2.fq.gz",
            },
        ]
        assert result == expected

    @patch("pathlib.Path.exists")
    def test_parse_fastq_filename_nonexistent_path(self, mock_exists):
        """测试不存在的路径"""
        mock_exists.return_value = False

        sample_path = Path("Sample-LIB001")
        result = self.processor.parse_fastq_filename(sample_path)

        assert result == []

    @patch("pathlib.Path.exists")
    @patch("pathlib.Path.glob")
    @patch("tqdm.tqdm")
    def test_build_libid_fastq_map(self, mock_tqdm, mock_glob, mock_exists):
        """测试构建libid映射"""
        mock_exists.return_value = True

        # 模拟Sample目录
        mock_sample_dir = Mock()
        mock_sample_dir.name = "Sample-LIB001"
        mock_glob.return_value = [mock_sample_dir]
        mock_tqdm.return_value = [mock_sample_dir]

        # 模拟parse_fastq_filename返回
        with patch.object(self.processor, "parse_fastq_filename") as mock_parse:
            mock_parse.return_value = [
                {"libid": "LIB001", "read_type": "R1", "path": "/path/R1.fastq.gz"}
            ]

            result = self.processor.build_libid_fastq_map(Path("test_path"))

            assert isinstance(result, pd.DataFrame)
            assert len(result) == 1
            assert result.iloc[0]["libid"] == "LIB001"

    @patch("pathlib.Path.exists")
    @patch("pandas.read_table")
    def test_read_or_build_config_existing_file(self, mock_read_table, mock_exists):
        """测试读取已存在的配置文件"""
        mock_exists.return_value = True

        # 模拟配置文件内容
        mock_df = pd.DataFrame(
            {
                "libid": ["LIB001"],
                "read_type": ["R1"],
                "path": ["/path/to/file.fastq.gz"],
            }
        )
        mock_read_table.return_value = mock_df

        result = self.processor.read_or_build_config(Path("test_dir"))

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 1
        mock_read_table.assert_called_once()

    @patch("pathlib.Path.exists")
    @patch("pandas.read_table")
    def test_read_or_build_config_rebuild_needed(self, mock_read_table, mock_exists):
        """测试需要重建配置文件的情况"""
        mock_exists.return_value = True

        # 模拟读取失败
        mock_read_table.side_effect = Exception("读取失败")

        with patch.object(self.processor, "build_libid_fastq_map") as mock_build:
            mock_build.return_value = pd.DataFrame(
                {
                    "libid": ["LIB001"],
                    "read_type": ["R1"],
                    "path": ["/path/to/file.fastq.gz"],
                }
            )

            result = self.processor.read_or_build_config(Path("test_dir"))

            assert isinstance(result, pd.DataFrame)
            mock_build.assert_called_once()

    @patch("pathlib.Path.glob")
    @patch("tqdm.tqdm")
    def test_load_config(self, mock_tqdm, mock_glob):
        """测试加载配置"""
        # 模拟目录结构
        mock_date_dir = Mock()
        mock_date_dir.is_dir.return_value = True
        mock_tcwl_dir = Mock()
        mock_tcwl_dir.name = "tcwl-test"
        mock_date_dir.glob.return_value = [mock_tcwl_dir]

        mock_glob.return_value = [mock_date_dir]  # base_dir.glob("20*")
        mock_tqdm.return_value = [mock_tcwl_dir]

        with patch.object(self.processor, "read_or_build_config") as mock_read:
            mock_read.return_value = pd.DataFrame(
                {
                    "libid": ["LIB001"],
                    "read_type": ["R1"],
                    "path": ["/path/to/file.fastq.gz"],
                }
            )

            result = self.processor.load_config(np.array(["tcwl-test"]))

            assert isinstance(result, pd.DataFrame)
            assert len(result) == 1


class TestScriptRunner:
    """测试ScriptRunner类"""

    def test_merge_or_link_command_single_file(self):
        """测试单文件命令生成"""
        fq_list = ["/path/to/file1.fastq.gz"]
        output_name = "output.fastq.gz"

        result = ScriptRunner.merge_or_link_command(fq_list, output_name)
        expected = "cp /path/to/file1.fastq.gz output.fastq.gz"

        assert result == expected

    def test_merge_or_link_command_multiple_files(self):
        """测试多文件命令生成"""
        fq_list = ["/path/to/file1.fastq.gz", "/path/to/file2.fastq.gz"]
        output_name = "output.fastq.gz"

        result = ScriptRunner.merge_or_link_command(fq_list, output_name)
        expected = (
            "cat /path/to/file1.fastq.gz /path/to/file2.fastq.gz > output.fastq.gz"
        )

        assert result == expected

    @patch("subprocess.run")
    def test_run_script_success(self, mock_run):
        """测试脚本执行成功"""
        mock_run.return_value = Mock(returncode=0)

        script_path = Path("test_script.sh")
        success, message = ScriptRunner.run_script(script_path)

        assert success is True
        assert "成功" in message
        mock_run.assert_called_once_with(
            ["bash", str(script_path)], check=True, capture_output=True, text=True
        )

    @patch("subprocess.run")
    def test_run_script_failure(self, mock_run):
        """测试脚本执行失败"""
        mock_run.side_effect = subprocess.CalledProcessError(
            1, "bash", stderr="Error message"
        )

        script_path = Path("test_script.sh")
        success, message = ScriptRunner.run_script(script_path)

        assert success is False
        assert "失败" in message
        assert "Error message" in message

    @patch("tqdm.tqdm")
    @patch("concurrent.futures.ThreadPoolExecutor")
    @patch("subprocess.run")  # 👈 最后 patch 的，最先传入
    @patch("pathlib.Path.glob")
    def test_run_scripts_in_parallel(
        self, mock_glob, mock_run, mock_executor, mock_tqdm
    ):
        """测试并行执行脚本"""
        # 模拟脚本列表
        script1 = Path("script1.sh")
        script2 = Path("script2.sh")
        mock_glob.return_value = [script1, script2]

        # mock run() 成功和失败
        def run_side_effect(args, check, capture_output, text):
            if "script1.sh" in args:
                return Mock(returncode=0, stdout="done")
            else:
                raise subprocess.CalledProcessError(1, args, stderr="Error")

        mock_run.side_effect = run_side_effect

        # ThreadPoolExecutor 模拟
        executor_instance = Mock()
        mock_executor.return_value.__enter__.return_value = executor_instance

        # 模拟 future 的 result
        future_success = Mock()
        future_success.result.return_value = (True, "成功: script1.sh")

        future_failed = Mock()
        future_failed.result.return_value = (False, "失败: script2.sh")

        # mock submit() 返回 future
        executor_instance.submit.side_effect = [future_success, future_failed]

        # mock as_completed
        with patch("concurrent.futures.as_completed") as mock_as_completed:
            mock_as_completed.return_value = [future_success, future_failed]
            mock_tqdm.return_value = [future_success, future_failed]

            result = ScriptRunner.run_scripts_in_parallel(Path("."))

        assert result["success"] == 1
        assert result["failed"] == 1


class TestWriteNextflowInput:
    """测试write_nextflow_input函数"""

    def setup_method(self):
        """每个测试方法前的设置"""
        self.temp_dir = Path(tempfile.mkdtemp())

    def teardown_method(self):
        """每个测试方法后的清理"""
        import shutil

        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_write_nextflow_input_empty_df(self):
        """测试空数据框"""
        empty_df = pd.DataFrame()

        result = write_nextflow_input(empty_df, self.temp_dir)

        assert result is None

    @patch("fq.project_fq.ScriptRunner.run_scripts_in_parallel")
    def test_write_nextflow_input_success(self, mock_run_scripts):
        """测试成功写入Nextflow输入"""
        # 准备测试数据
        test_df = pd.DataFrame(
            {
                "sample_id": ["SAMPLE1", "SAMPLE1", "SAMPLE2"],
                "read_type": ["R1", "R2", "R1"],
                "path": [
                    "/path/file1_R1.fq.gz",
                    "/path/file1_R2.fq.gz",
                    "/path/file2_R1.fq.gz",
                ],
            }
        )

        mock_run_scripts.return_value = {"success": 2, "failed": 0}

        result = write_nextflow_input(test_df, self.temp_dir)

        assert result is not None
        assert result["success"] == 2
        assert result["failed"] == 0

        # 检查脚本目录是否创建
        scripts_dir = self.temp_dir / "scripts"
        assert scripts_dir.exists()

        mock_run_scripts.assert_called_once()

    def test_write_nextflow_input_with_na_paths(self):
        """测试包含NA路径的数据"""
        test_df = pd.DataFrame(
            {
                "sample_id": ["SAMPLE1", "SAMPLE1"],
                "read_type": ["R1", "R2"],
                "path": ["/path/file1_R1.fq.gz", None],
            }
        )

        result = write_nextflow_input(test_df, self.temp_dir)

        # 应该只生成一个脚本文件
        scripts_dir = self.temp_dir / "scripts"
        if scripts_dir.exists():
            scripts = list(scripts_dir.glob("*.sh"))
            assert len(scripts) <= 1


class TestValidateSampleInfo:
    """测试validate_sample_info函数"""

    def test_validate_sample_info_valid(self):
        """测试有效的样品信息"""
        valid_df = pd.DataFrame(
            {
                "libid": ["LIB001", "LIB002"],
                "sample_id": ["SAMPLE1", "SAMPLE2"],
                "dir_name": ["dir1", "dir2"],
                "extra_col": ["extra1", "extra2"],
            }
        )

        result = validate_sample_info(valid_df)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2

    def test_validate_sample_info_missing_columns(self):
        """测试缺少必要列"""
        invalid_df = pd.DataFrame(
            {
                "libid": ["LIB001"],
                "sample_id": ["SAMPLE1"],
                # 缺少 dir_name 列
            }
        )

        with pytest.raises(ValueError, match="样品信息文件缺少必要列"):
            validate_sample_info(invalid_df)

    def test_validate_sample_info_duplicates(self):
        """测试去重功能"""
        df_with_duplicates = pd.DataFrame(
            {
                "libid": ["LIB001", "LIB001", "LIB002"],
                "sample_id": ["SAMPLE1", "SAMPLE1", "SAMPLE2"],
                "dir_name": ["dir1", "dir1", "dir2"],
            }
        )

        result = validate_sample_info(df_with_duplicates)

        assert len(result) == 2  # 应该去除一个重复行

    def test_validate_sample_info_with_nulls(self):
        """测试包含空值的数据"""
        df_with_nulls = pd.DataFrame(
            {
                "libid": ["LIB001", None],
                "sample_id": ["SAMPLE1", "SAMPLE2"],
                "dir_name": ["dir1", "dir2"],
            }
        )

        # 应该能够处理空值而不抛出异常
        result = validate_sample_info(df_with_nulls)
        assert isinstance(result, pd.DataFrame)


class TestLogStatistics:
    """测试log_statistics函数"""

    @patch("fq.project_fq.logger")
    def test_log_statistics_empty_df(self, mock_logger):
        """测试空数据框的统计"""
        empty_df = pd.DataFrame()

        log_statistics(empty_df)

        mock_logger.error.assert_called_with("数据框为空")

    @patch("fq.project_fq.logger")
    def test_log_statistics_complete_data(self, mock_logger):
        """测试完整数据的统计"""
        complete_df = pd.DataFrame(
            {
                "sample_id": ["SAMPLE1", "SAMPLE1", "SAMPLE2"],
                "libid": ["LIB001", "LIB002", "LIB003"],
                "path": ["/path1", "/path2", "/path3"],
            }
        )

        log_statistics(complete_df)

        # 检查是否记录了成功信息
        mock_logger.success.assert_called_with("所有样品数据已找到")

    @patch("fq.project_fq.logger")
    def test_log_statistics_missing_data(self, mock_logger):
        """测试缺失数据的统计"""
        missing_df = pd.DataFrame(
            {
                "sample_id": ["SAMPLE1", "SAMPLE2"],
                "libid": ["LIB001", "LIB002"],
                "path": ["/path1", None],
            }
        )

        log_statistics(missing_df)

        # 检查是否记录了错误信息
        mock_logger.error.assert_called()


class TestConstants:
    """测试常量定义"""

    def test_fastq_extensions(self):
        """测试FASTQ扩展名常量"""
        assert "*.fastq.gz" in FASTQ_EXTENSIONS
        assert "*.fq.gz" in FASTQ_EXTENSIONS
        assert len(FASTQ_EXTENSIONS) == 2

    def test_read_type_patterns(self):
        """测试读取类型模式常量"""
        assert READ_TYPE_PATTERNS["combined_R1.fastq.gz"] == "R1"
        assert READ_TYPE_PATTERNS["combined_R2.fastq.gz"] == "R2"
        assert READ_TYPE_PATTERNS["_R1.fastq.gz"] == "R1"
        assert READ_TYPE_PATTERNS["_R2.fastq.gz"] == "R2"
        assert READ_TYPE_PATTERNS["_1.fastq.gz"] == "R1"
        assert READ_TYPE_PATTERNS["_2.fastq.gz"] == "R2"


# 集成测试
class TestIntegration:
    """集成测试"""

    def setup_method(self):
        """每个测试方法前的设置"""
        self.temp_dir = Path(tempfile.mkdtemp())
        self.processor = FastqProcessor(self.temp_dir)

    def teardown_method(self):
        """每个测试方法后的清理"""
        import shutil

        shutil.rmtree(self.temp_dir, ignore_errors=True)

    @patch("pathlib.Path.glob")
    @patch("pandas.DataFrame.to_csv")
    def test_full_workflow_mock(self, mock_to_csv, mock_glob):
        """测试完整工作流程（使用mock）"""
        # 模拟样品信息
        sample_df = pd.DataFrame(
            {
                "libid": ["LIB001", "LIB002"],
                "sample_id": ["SAMPLE1", "SAMPLE2"],
                "dir_name": ["tcwl-test1", "tcwl-test2"],
            }
        )

        # 模拟目录结构
        mock_date_dir = Mock()
        mock_date_dir.is_dir.return_value = True
        mock_tcwl_dir = Mock()
        mock_tcwl_dir.name = "tcwl-test1"
        mock_date_dir.glob.return_value = [mock_tcwl_dir]
        mock_glob.return_value = [mock_date_dir]

        # 模拟配置加载
        with patch.object(self.processor, "read_or_build_config") as mock_read_config:
            mock_read_config.return_value = pd.DataFrame(
                {
                    "libid": ["LIB001"],
                    "read_type": ["R1"],
                    "path": ["/path/to/file.fastq.gz"],
                    "dir_name": ["tcwl-test1"],
                }
            )

            # 执行配置加载
            result = self.processor.load_config(sample_df["dir_name"].unique())

            # 验证结果
            assert isinstance(result, pd.DataFrame)
            assert len(result) > 0

            # 验证合并操作
            merged_df = sample_df.merge(result, on="libid", how="left")
            assert "path" in merged_df.columns


# 性能测试
class TestPerformance:
    """性能测试"""

    def test_large_dataframe_processing(self):
        """测试大数据框处理性能"""
        # 创建大数据框
        large_df = pd.DataFrame(
            {
                "libid": [f"LIB{i:06d}" for i in range(1000)],
                "sample_id": [f"SAMPLE{i:06d}" for i in range(1000)],
                "dir_name": [f"dir{i}" for i in range(1000)],
            }
        )

        # 测试验证函数的性能
        import time

        start_time = time.time()

        result = validate_sample_info(large_df)

        end_time = time.time()
        processing_time = end_time - start_time

        # 验证结果
        assert len(result) == 1000
        assert processing_time < 1.0  # 应该在1秒内完成


# 参数化测试
class TestParametrized:
    """参数化测试"""

    @pytest.mark.parametrize(
        "filename,expected",
        [
            ("sample_combined_R1.fastq.gz", "R1"),
            ("sample_combined_R2.fastq.gz", "R2"),
            ("sample_R1.fastq.gz", "R1"),
            ("sample_R2.fastq.gz", "R2"),
            ("sample_1.fastq.gz", "R1"),
            ("sample_2.fastq.gz", "R2"),
            ("sample_unknown.fastq.gz", None),
            ("sample.txt", None),
        ],
    )
    def test_determine_read_type_parametrized(self, filename, expected):
        """参数化测试读取类型识别"""
        temp_dir = Path(tempfile.mkdtemp())
        processor = FastqProcessor(temp_dir)

        result = processor._determine_read_type(filename)
        assert result == expected

        # 清理
        import shutil

        shutil.rmtree(temp_dir, ignore_errors=True)

    @pytest.mark.parametrize(
        "input_files,expected_command_type",
        [
            (["/path/file1.fastq.gz"], "cp"),
            (["/path/file1.fastq.gz", "/path/file2.fastq.gz"], "cat"),
            (["/path/f1.fastq.gz", "/path/f2.fastq.gz", "/path/f3.fastq.gz"], "cat"),
        ],
    )
    def test_merge_command_parametrized(self, input_files, expected_command_type):
        """参数化测试合并命令生成"""
        output_name = "output.fastq.gz"

        result = ScriptRunner.merge_or_link_command(input_files, output_name)

        assert result.startswith(expected_command_type)
        assert output_name in result


# 错误处理测试
class TestErrorHandling:
    """错误处理测试"""

    def test_file_permission_error(self):
        """测试文件权限错误处理"""
        temp_dir = Path(tempfile.mkdtemp())
        processor = FastqProcessor(temp_dir)

        # 模拟权限错误
        with patch("pathlib.Path.glob") as mock_glob:
            mock_glob.side_effect = PermissionError("Permission denied")

            result = processor.build_libid_fastq_map(temp_dir)

            # 应该返回空的DataFrame而不是抛出异常
            assert isinstance(result, pd.DataFrame)
            assert len(result) == 0

        # 清理
        import shutil

        shutil.rmtree(temp_dir, ignore_errors=True)

    def test_malformed_data_handling(self):
        """测试格式错误的数据处理"""
        malformed_df = pd.DataFrame(
            {
                "libid": ["LIB001", "LIB002", ""],
                "sample_id": ["SAMPLE1", "", "SAMPLE3"],
                "dir_name": ["dir1", "dir2", "dir3"],
            }
        )

        # 应该能够处理而不抛出异常
        result = validate_sample_info(malformed_df)
        assert isinstance(result, pd.DataFrame)


if __name__ == "__main__":
    # 运行测试的示例命令
    # pytest test_fastq_processor.py -v
    # pytest test_fastq_processor.py::TestFastqProcessor::test_init_valid_directory -v
    # pytest test_fastq_processor.py -k "test_determine_read_type" -v
    # pytest test_fastq_processor.py --cov=fastq_processor --cov-report=html
    pass
