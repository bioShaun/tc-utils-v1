from pathlib import Path

import pandas as pd
import pytest

from chip.panel2bed import load_split_bed, main, split_bed_dataframe


def write_basic_inputs(tmp_path: Path) -> tuple[Path, Path]:
    """写入最小输入文件：design.tsv 与 genome.fai。"""
    design_table = tmp_path / "design.tsv"
    genome_fai = tmp_path / "genome.fa.fai"

    design_df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "pos": [100, 400],
            "probe_start": [80, 380],
            "probe_end": [120, 420],
        }
    )
    design_df.to_csv(design_table, sep="\t", index=False)
    genome_fai.write_text("chr1\t1000\t0\t0\t0\n", encoding="utf-8")
    return design_table, genome_fai


def read_lines(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8").strip().splitlines()


def test_panel2bed_without_split_bed_outputs_original_files(tmp_path: Path) -> None:
    design_table, genome_fai = write_basic_inputs(tmp_path)
    out_dir = tmp_path / "out"

    main(
        design_table=design_table,
        genome_fai=genome_fai,
        probe_id="panel_v1",
        out_path=out_dir,
        flank_size=200,
        split_bed=None,
    )

    assert (out_dir / "panel_v1.id").exists()
    assert (out_dir / "panel_v1.bed").exists()
    assert (out_dir / "panel_v1.snpcalling.bed").exists()
    assert not (out_dir / "panel_v1.split.bed").exists()
    assert not (out_dir / "panel_v1.snpcalling.split.bed").exists()

    assert read_lines(out_dir / "panel_v1.bed") == [
        "chr1\t99\t100",
        "chr1\t399\t400",
    ]
    assert read_lines(out_dir / "panel_v1.snpcalling.bed") == [
        "chr1\t0\t200",
        "chr1\t300\t500",
    ]


def test_panel2bed_with_split_bed_outputs_split_files_only(tmp_path: Path) -> None:
    design_table, genome_fai = write_basic_inputs(tmp_path)
    split_bed = tmp_path / "split.bed"
    out_dir = tmp_path / "out"
    split_bed.write_text(
        "chr1\t0\t300\tchr1_part1\n"
        "chr1\t300\t600\tchr1_part2\n",
        encoding="utf-8",
    )
    # 不显式传 --split-genome-fai，使用同目录默认文件
    (tmp_path / "split.genome.fa.fai").write_text(
        "chr1_part1\t300\t0\t0\t0\n"
        "chr1_part2\t300\t0\t0\t0\n",
        encoding="utf-8",
    )

    main(
        design_table=design_table,
        genome_fai=genome_fai,
        probe_id="panel_v1",
        out_path=out_dir,
        flank_size=200,
        split_bed=split_bed,
    )

    assert (out_dir / "panel_v1.id").exists()
    assert (out_dir / "panel_v1.split.bed").exists()
    assert (out_dir / "panel_v1.snpcalling.split.bed").exists()
    assert not (out_dir / "panel_v1.bed").exists()
    assert not (out_dir / "panel_v1.snpcalling.bed").exists()

    assert read_lines(out_dir / "panel_v1.split.bed") == [
        "chr1_part1\t99\t100",
        "chr1_part2\t99\t100",
    ]
    assert read_lines(out_dir / "panel_v1.snpcalling.split.bed") == [
        "chr1_part1\t0\t200",
        "chr1_part2\t0\t200",
    ]


def test_split_logic_matches_reference_filter_rule() -> None:
    bed_df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1"],
            "start": [0, 49, 50],
            "end": [10, 60, 70],
        }
    )
    split_df = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "split_start": [0],
            "split_end": [50],
            "new_chrom": ["chr1_part1"],
        }
    )

    result_df = split_bed_dataframe(
        bed_df=bed_df,
        split_bed_df=split_df,
        start_col="start",
        end_col="end",
    )

    assert result_df.to_dict(orient="records") == [
        {"new_chrom": "chr1_part1", "new_start": 0, "new_end": 10},
        {"new_chrom": "chr1_part1", "new_start": 49, "new_end": 60},
    ]


def test_split_bed_missing_or_invalid_raises_error(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        load_split_bed(tmp_path / "missing_split.bed")

    invalid_split_bed = tmp_path / "invalid_split.bed"
    invalid_split_bed.write_text("chr1\t0\t100\n", encoding="utf-8")
    with pytest.raises(ValueError):
        load_split_bed(invalid_split_bed)


def test_split_output_sorted_by_split_genome_fai(tmp_path: Path) -> None:
    design_table = tmp_path / "design.tsv"
    design_df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "pos": [100, 100],
            "probe_start": [80, 80],
            "probe_end": [120, 120],
        }
    )
    design_df.to_csv(design_table, sep="\t", index=False)

    genome_fai = tmp_path / "genome.fa.fai"
    genome_fai.write_text(
        "chr1\t1000\t0\t0\t0\n"
        "chr2\t1000\t0\t0\t0\n",
        encoding="utf-8",
    )

    split_bed = tmp_path / "split.bed"
    split_bed.write_text(
        "chr1\t0\t300\tchr1_part2\n"
        "chr2\t0\t300\tchr2_part1\n",
        encoding="utf-8",
    )
    split_genome_fai = tmp_path / "custom.split.genome.fa.fai"
    split_genome_fai.write_text(
        "chr2_part1\t300\t0\t0\t0\n"
        "chr1_part2\t300\t0\t0\t0\n",
        encoding="utf-8",
    )

    out_dir = tmp_path / "out"
    main(
        design_table=design_table,
        genome_fai=genome_fai,
        probe_id="panel_v1",
        out_path=out_dir,
        flank_size=200,
        split_bed=split_bed,
        split_genome_fai=split_genome_fai,
    )

    assert read_lines(out_dir / "panel_v1.split.bed") == [
        "chr2_part1\t99\t100",
        "chr1_part2\t99\t100",
    ]


def test_split_bed_without_split_genome_fai_raises_error(tmp_path: Path) -> None:
    design_table, genome_fai = write_basic_inputs(tmp_path)
    split_bed = tmp_path / "split.bed"
    split_bed.write_text(
        "chr1\t0\t300\tchr1_part1\n"
        "chr1\t300\t600\tchr1_part2\n",
        encoding="utf-8",
    )
    out_dir = tmp_path / "out"

    with pytest.raises(ValueError, match="split.genome.fa.fai"):
        main(
            design_table=design_table,
            genome_fai=genome_fai,
            probe_id="panel_v1",
            out_path=out_dir,
            flank_size=200,
            split_bed=split_bed,
            split_genome_fai=None,
        )
