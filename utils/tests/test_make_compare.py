from pathlib import Path

import pytest
from typer.testing import CliRunner

from utils.make_compare import app, generate_pairs, read_sample_list

runner = CliRunner()


def test_read_sample_list(tmp_path):
    # Create a temporary sample file
    sample_file = tmp_path / "samples.txt"
    sample_file.write_text("sample1\nsample2\nsample3\n")

    samples = read_sample_list(sample_file)
    assert len(samples) == 3
    assert samples == ["sample1", "sample2", "sample3"]


def test_read_sample_list_empty_lines(tmp_path):
    # Test with empty lines and whitespace
    sample_file = tmp_path / "samples.txt"
    sample_file.write_text("sample1\n\n  sample2  \n\nsample3\n")

    samples = read_sample_list(sample_file)
    assert len(samples) == 3
    assert samples == ["sample1", "sample2", "sample3"]


def test_generate_pairs():
    samples = ["A", "B", "C"]
    pairs = generate_pairs(samples)

    assert len(pairs) == 3  # 3 choose 2 = 3 combinations
    assert ("A", "B") in pairs
    assert ("A", "C") in pairs
    assert ("B", "C") in pairs


def test_generate_pairs_two_samples():
    samples = ["A", "B"]
    pairs = generate_pairs(samples)

    assert len(pairs) == 1
    assert pairs == [("A", "B")]


def test_generate_pairs_one_sample():
    samples = ["A"]
    pairs = generate_pairs(samples)

    assert len(pairs) == 0  # Can't generate pairs with only one sample


def test_cli_basic_usage(tmp_path):
    # Create a temporary sample file
    sample_file = tmp_path / "samples.txt"
    sample_file.write_text("A\nB\nC\n")

    result = runner.invoke(app, [str(sample_file)])
    assert result.exit_code == 0
    assert "Found 3 samples, generating 3 pairwise combinations" in result.stdout
    assert "A\tB" in result.stdout
    assert "A\tC" in result.stdout
    assert "B\tC" in result.stdout


def test_cli_output_file(tmp_path):
    # Create input and output files
    input_file = tmp_path / "samples.txt"
    output_file = tmp_path / "results.txt"
    input_file.write_text("A\nB\nC\n")

    result = runner.invoke(app, ["--output", str(output_file), str(input_file)])
    assert result.exit_code == 0
    assert f"Results saved to {output_file}" in result.stdout

    # Check output file contents
    content = output_file.read_text()
    assert "Found 3 samples, generating 3 pairwise combinations" in content
    assert "A\tB" in content
    assert "A\tC" in content
    assert "B\tC" in content


def test_cli_empty_file(tmp_path):
    # Create an empty file
    sample_file = tmp_path / "empty.txt"
    sample_file.write_text("")

    result = runner.invoke(app, [str(sample_file)])
    assert result.exit_code == 1
    assert "No samples found in the input file" in result.stderr


def test_cli_nonexistent_file():
    result = runner.invoke(app, ["nonexistent.txt"])
    assert result.exit_code == 1
    assert "does not exist" in result.stderr


def test_cli_help():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "Generate all possible pairwise combinations" in result.stdout
    assert "INPUT_FILE" in result.stdout
    assert "--output" in result.stdout
