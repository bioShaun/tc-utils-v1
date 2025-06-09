import gzip
import logging
from pathlib import Path

import polars as pl
import pytest

from panel.cdsCovEvaluation_bamdst_02x_polars import (
    DEFAULT_DEPTH_THRESHOLD,
    DEFAULT_FLOAT_PRECISION,
    CoverageAnalysisError,
    DataValidationError,
    FileProcessingError,
    calculate_coverage_ratio,
    load_bed_files,
    load_sample_list,
    load_single_depth_file,
    main,
    merge_chr,
    validate_dataframe,
    validate_input_path,
    write_output,
)


# ---- 日志桥接 fixture ----
@pytest.fixture(autouse=True)
def caplog_bridge(caplog):
    logging.getLogger().setLevel(logging.DEBUG)
    yield


# ========== Exception Tests ==========
class TestCustomExceptions:
    def test_coverage_analysis_error(self):
        with pytest.raises(CoverageAnalysisError):
            raise CoverageAnalysisError("Test error")

    def test_data_validation_error(self):
        with pytest.raises(DataValidationError):
            raise DataValidationError("Validation error")

    def test_file_processing_error(self):
        with pytest.raises(FileProcessingError):
            raise FileProcessingError("File processing error")


# ========== Path Validation ==========
class TestValidateInputPath:
    def test_validate_existing_directory(self, tmp_path):
        validate_input_path(tmp_path, "directory")

    def test_validate_existing_file(self, tmp_path):
        f = tmp_path / "a.txt"
        f.write_text("abc")
        validate_input_path(f, "file")

    def test_validate_nonexistent_path(self, tmp_path):
        with pytest.raises(DataValidationError):
            validate_input_path(tmp_path / "notexists", "file")

    def test_validate_file_as_directory(self, tmp_path):
        f = tmp_path / "a.txt"
        f.write_text("abc")
        with pytest.raises(DataValidationError):
            validate_input_path(f, "directory")

    def test_validate_directory_as_file(self, tmp_path):
        d = tmp_path / "d"
        d.mkdir()
        with pytest.raises(DataValidationError):
            validate_input_path(d, "file")


# ========== DataFrame Validation ==========
class TestValidateDataframe:
    def test_validate_empty_dataframe(self, caplog):
        empty_df = pl.DataFrame()
        validate_dataframe(empty_df, "test_df")
        assert "test_df is empty" in caplog.text

    def test_validate_dataframe_with_required_columns(self):
        df = pl.DataFrame({"col1": [1, 2], "col2": [3, 4]})
        validate_dataframe(df, "test_df", ["col1", "col2"])

    def test_validate_dataframe_missing_columns(self):
        df = pl.DataFrame({"col1": [1, 2]})
        with pytest.raises(DataValidationError):
            validate_dataframe(df, "test_df", ["col1", "col2"])

    def test_validate_nonempty_dataframe_no_required_columns(self):
        df = pl.DataFrame({"col1": [1, 2]})
        validate_dataframe(df, "test_df")


# ========== merge_chr ==========
class TestMergeChr:
    def create_test_split_bed(self, tmp_path):
        p = tmp_path / "split.bed"
        p.write_text("chr1_split\t0\t1000\tchr1\nchr2_split\t0\t2000\tchr2\n")
        return p

    def test_merge_chr_success(self, tmp_path):
        df = pl.DataFrame(
            {"chrom": ["chr1", "chr2"], "start": [100, 200], "end": [200, 300]}
        )
        split_bed_file = self.create_test_split_bed(tmp_path)
        result = merge_chr(df, split_bed_file)
        assert result.shape[0] == 2
        assert (
            "chrom" in result.columns
            and "start" in result.columns
            and "end" in result.columns
        )
        assert result["start"].to_list() == [100, 200]
        assert result["end"].to_list() == [200, 300]

    def test_merge_chr_no_matches(self, tmp_path, caplog):
        df = pl.DataFrame({"chrom": ["chr3"], "start": [1], "end": [2]})
        split_bed_file = self.create_test_split_bed(tmp_path)
        result = merge_chr(df, split_bed_file)
        assert result.equals(df)
        assert "No matching chromosomes found" in caplog.text

    def test_merge_chr_invalid_file(self, tmp_path):
        df = pl.DataFrame({"chrom": ["chr1"], "start": [1], "end": [2]})
        with pytest.raises(DataValidationError):
            merge_chr(df, tmp_path / "notexists.bed")


# ========== load_single_depth_file ==========
class TestLoadSingleDepthFile:
    def create_depth_file(self, tmp_path, lines):
        import gzip

        f = tmp_path / "depth.tsv.gz"
        with gzip.open(f, "wt") as fp:
            fp.write("\n".join(lines))
        return f

    def test_load_single_depth_file_success(self, tmp_path):
        content = ["chr1\t100\t10", "chr1\t101\t20", "chr1\t102\t5"]
        f = self.create_depth_file(tmp_path, content)
        result = load_single_depth_file(f, "sample1")
        assert result is not None
        assert "sample1" in result.columns and result.shape[0] == 3
        mean = sum([10, 20, 5]) / 3
        threshold = mean * DEFAULT_DEPTH_THRESHOLD
        expected = [d >= threshold for d in [10, 20, 5]]
        assert result["sample1"].to_list() == expected

    def test_load_single_depth_file_empty(self, tmp_path, caplog):
        f = self.create_depth_file(tmp_path, [])
        result = load_single_depth_file(f, "sample1")
        assert result is None
        assert "Empty depth file" in caplog.text

    def test_load_single_depth_file_zero_depth(self, tmp_path, caplog):
        f = self.create_depth_file(tmp_path, ["chr1\t100\t0", "chr1\t101\t0"])
        result = load_single_depth_file(f, "sample1")
        assert result is not None
        assert "Zero or null mean depth" in caplog.text
        # 所有值都应该是 True（因为深度阈值为0，而深度也为0）
        assert all(x is True for x in result["sample1"].to_list())


# ========== load_bed_files ==========
class TestLoadBedFiles:
    def create_bed_structure(self, tmp_path, samples=("s1", "s2"), nrows=3):
        for s in samples:
            d = tmp_path / s
            d.mkdir()
            lines = [f"chr1\t{100+i}\t{10*(i+1)}" for i in range(nrows)]
            with gzip.open(d / "depth.tsv.gz", "wt") as fp:
                fp.write("\n".join(lines))
        return tmp_path

    def test_load_bed_files_success(self, tmp_path):
        bed_dir = self.create_bed_structure(tmp_path)
        bed_df, matrix_df = load_bed_files(bed_dir)
        assert not bed_df.is_empty()
        assert not matrix_df.is_empty()
        assert bed_df.shape[0] == 3
        assert matrix_df.shape[1] == 2
        assert set(bed_df.columns) == {"chrom", "start", "end"}

    def test_load_bed_files_with_sample_list(self, tmp_path):
        bed_dir = self.create_bed_structure(tmp_path, samples=("a", "b", "c"))
        bed_df, matrix_df = load_bed_files(bed_dir, sample_list=["a"])
        assert not bed_df.is_empty()
        assert matrix_df.shape[1] == 1

    def test_load_bed_files_no_files(self, tmp_path):
        empty_dir = tmp_path / "e"
        empty_dir.mkdir()
        with pytest.raises(FileProcessingError):
            load_bed_files(empty_dir)

    def test_load_bed_files_invalid_directory(self, tmp_path):
        with pytest.raises(DataValidationError):
            load_bed_files(tmp_path / "notexists")


# ========== calculate_coverage_ratio ==========
class TestCalculateCoverageRatio:
    def test_calculate_coverage_ratio_success(self):
        df = pl.DataFrame({"a": [True, False, True], "b": [True, True, False]})
        result = calculate_coverage_ratio(df)
        assert result.to_list() == [1.0, 0.5, 0.5]

    def test_calculate_coverage_ratio_empty(self):
        df = pl.DataFrame()
        result = calculate_coverage_ratio(df)
        assert result.to_list() == []

    def test_calculate_coverage_ratio_single_sample(self):
        df = pl.DataFrame({"a": [True, False, True]})
        assert calculate_coverage_ratio(df).to_list() == [1.0, 0.0, 1.0]


# ========== load_sample_list ==========
class TestLoadSampleList:
    def test_load_sample_list_success(self, tmp_path):
        p = tmp_path / "samples.txt"
        p.write_text("a\nb\nc\n")
        assert load_sample_list(p) == ["a", "b", "c"]

    def test_load_sample_list_empty_file(self, tmp_path):
        p = tmp_path / "empty.txt"
        p.write_text("")
        assert load_sample_list(p) == []

    def test_load_sample_list_nonexistent_file(self, tmp_path):
        with pytest.raises(FileProcessingError):
            load_sample_list(tmp_path / "notexists.txt")


# ========== write_output ==========
class TestWriteOutput:
    def test_write_output_success(self, tmp_path):
        df = pl.DataFrame(
            {
                "chrom": ["chr1", "chr2"],
                "start": [100, 200],
                "end": [200, 300],
                "coverage_0.2x": [0.75, 0.5],
            }
        )
        out = tmp_path / "outdir" / "res.tsv"
        write_output(df, out)
        assert out.exists()
        df2 = pl.read_csv(out, separator="\t")
        assert df2.shape == df.shape and set(df2.columns) == set(df.columns)

    def test_write_output_creates_directory(self, tmp_path):
        df = pl.DataFrame({"col1": [1, 2]})
        out = tmp_path / "a" / "b" / "f.tsv"
        write_output(df, out)
        assert out.exists() and out.parent.exists()


# ========== Main Typer CLI 测试 ==========
import typer
from typer.testing import CliRunner


class TestMain:
    def create_complete_test_structure(self, tmp_path):
        cds_dir = tmp_path / "cds_cov"
        for s in ["s1", "s2"]:
            d = cds_dir / s
            d.mkdir(parents=True)
            lines = ["chr1\t100\t10", "chr1\t200\t20", "chr1\t300\t15"]
            with gzip.open(d / "depth.tsv.gz", "wt") as fp:
                fp.write("\n".join(lines))
        sample_file = tmp_path / "samples.txt"
        sample_file.write_text("s1\ns2\n")
        split_bed_file = tmp_path / "split.bed"
        split_bed_file.write_text("chr1_new\t0\t1000\tchr1\n")
        return cds_dir, sample_file, split_bed_file

    def test_main_basic_functionality(self, tmp_path):
        cds_dir, sample_file, split_bed_file = self.create_complete_test_structure(
            tmp_path
        )
        output_file = tmp_path / "out.tsv"
        app = typer.Typer()
        app.command()(main)
        runner = CliRunner()
        result = runner.invoke(
            app, [str(cds_dir), str(output_file), "--log-level", "ERROR"]
        )
        assert result.exit_code == 0 and output_file.exists()

    def test_main_with_sample_list(self, tmp_path):
        cds_dir, sample_file, split_bed_file = self.create_complete_test_structure(
            tmp_path
        )
        output_file = tmp_path / "out.tsv"
        app = typer.Typer()
        app.command()(main)
        runner = CliRunner()
        result = runner.invoke(
            app,
            [
                str(cds_dir),
                str(output_file),
                "--sample-path",
                str(sample_file),
                "--log-level",
                "ERROR",
            ],
        )
        assert result.exit_code == 0 and output_file.exists()

    def test_main_with_split_bed(self, tmp_path):
        cds_dir, sample_file, split_bed_file = self.create_complete_test_structure(
            tmp_path
        )
        output_file = tmp_path / "out.tsv"
        app = typer.Typer()
        app.command()(main)
        runner = CliRunner()
        result = runner.invoke(
            app,
            [
                str(cds_dir),
                str(output_file),
                "--split-bed",
                str(split_bed_file),
                "--log-level",
                "ERROR",
            ],
        )
        assert result.exit_code == 0 and output_file.exists()

    def test_main_invalid_input_directory(self, tmp_path):
        output_file = tmp_path / "o.tsv"
        app = typer.Typer()
        app.command()(main)
        runner = CliRunner()
        result = runner.invoke(
            app, [str(tmp_path / "notexists"), str(output_file), "--log-level", "ERROR"]
        )
        assert result.exit_code == 1

    def test_main_no_depth_files(self, tmp_path):
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        output_file = tmp_path / "out.tsv"
        app = typer.Typer()
        app.command()(main)
        runner = CliRunner()
        result = runner.invoke(
            app, [str(empty_dir), str(output_file), "--log-level", "ERROR"]
        )
        assert result.exit_code == 1


# ========== 常量 ==========
class TestConstants:

    def test_default_depth_threshold(self):
        assert DEFAULT_DEPTH_THRESHOLD == 0.2

    def test_default_float_precision(self):
        assert DEFAULT_FLOAT_PRECISION == 3


# ========== 集成测试 ==========
class TestIntegration:
    def test_full_workflow(self, tmp_path):
        cds_dir = tmp_path / "cov"
        for s in ["a", "b"]:
            d = cds_dir / s
            d.mkdir(parents=True)
            with gzip.open(d / "depth.tsv.gz", "wt") as fp:
                fp.write("\n".join([f"chr1\t{i*10}\t{10+i}" for i in range(5)]))

        sample_file = tmp_path / "samples.txt"
        sample_file.write_text("a\nb\n")
        split_bed_file = tmp_path / "split.bed"
        split_bed_file.write_text("chr1_merge\t0\t1000\tchr1\n")
        out = tmp_path / "final.tsv"
        app = typer.Typer()
        app.command()(main)
        runner = CliRunner()
        result = runner.invoke(
            app,
            [
                str(cds_dir),
                str(out),
                "--sample-path",
                str(sample_file),
                "--split-bed",
                str(split_bed_file),
                "--log-level",
                "INFO",
            ],
        )
        assert result.exit_code == 0 and out.exists()
        df = pl.read_csv(out, separator="\t")
        assert "chrom" in df.columns and "coverage_0.2x" in df.columns


# ========== 性能测试（可选，慢用slow标记） ==========
import time


@pytest.mark.slow
def test_large_dataset_performance(tmp_path):
    cds_dir = tmp_path / "large"
    num_samples = 8
    num_regions = 300
    for i in range(num_samples):
        d = cds_dir / f"s{i}"
        d.mkdir(parents=True)
        content = [f"chr1\t{j*10}\t{10 + (j+i)%30}" for j in range(num_regions)]
        with gzip.open(d / "depth.tsv.gz", "wt") as fp:
            fp.write("\n".join(content))
    out = tmp_path / "o.tsv"
    app = typer.Typer()
    app.command()(main)
    runner = CliRunner()
    t0 = time.time()
    result = runner.invoke(app, [str(cds_dir), str(out), "--log-level", "ERROR"])
    t1 = time.time()
    assert result.exit_code == 0 and out.exists() and t1 - t0 < 30


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
