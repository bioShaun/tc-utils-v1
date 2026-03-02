from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas as pd
import pytest
from pydantic import ValidationError
from typer.testing import CliRunner

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "probe-va-overlap-filter-ai.py"
MODULE_SPEC = spec_from_file_location("probe_va_overlap_filter_ai", SCRIPT_PATH)
if MODULE_SPEC is None or MODULE_SPEC.loader is None:
    raise ImportError(f"无法加载脚本模块: {SCRIPT_PATH}")

probe_module = module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(probe_module)

runner = CliRunner()


def test_processing_config_rejects_invalid_values() -> None:
    with pytest.raises(ValidationError):
        probe_module.ProcessingConfig(threads=0, variant_cutoff=3, indel_cutoff=0)

    with pytest.raises(ValidationError):
        probe_module.ProcessingConfig(threads=1, variant_cutoff=-1, indel_cutoff=0)


def test_write_filtered_output_in_chunks_filters_rows(tmp_path: Path) -> None:
    ann_table = tmp_path / "ann.tsv"
    ann_table.write_text(
        "chrom\tprobe_start\tprobe_end\tid\tgene\n"
        "chr1\t100\t120\tp1\tG1\n"
        "chr1\t130\t150\tp2\tG2\n"
        "chr1\t160\t180\tp3\tG3\n",
        encoding="utf-8",
    )
    out_table = tmp_path / "out.tsv"
    id_list = tmp_path / "ids.txt"
    id_list.write_text("p1\np3\n", encoding="utf-8")

    ann_processor = probe_module.AnnDataProcessor(ann_table=ann_table)
    config = probe_module.ProcessingConfig(threads=2, variant_cutoff=3, indel_cutoff=1)
    variant_overlap_df = pd.DataFrame(
        {"id": ["p1", "p2", "p3"], "variant_overlap": [2, 5, 1]}
    )
    indel_overlap_df = pd.DataFrame(
        {"id": ["p1", "p2", "p3"], "indel_overlap": [0, 0, 2]}
    )

    probe_module.write_filtered_output_in_chunks(
        ann_data_processor=ann_processor,
        out_table=out_table,
        config=config,
        variant_overlap_df=variant_overlap_df,
        indel_overlap_df=indel_overlap_df,
        id_list=id_list,
    )

    result_df = pd.read_table(out_table)
    assert list(result_df.columns) == [
        "chrom",
        "probe_start",
        "probe_end",
        "id",
        "gene",
        "indel_overlap",
        "variant_overlap",
    ]
    assert result_df["id"].tolist() == ["p1"]
    assert result_df["variant_overlap"].tolist() == [2]
    assert result_df["indel_overlap"].tolist() == [0]


def test_write_filtered_output_in_chunks_writes_header_when_empty(tmp_path: Path) -> None:
    ann_table = tmp_path / "ann.tsv"
    ann_table.write_text(
        "chrom\tprobe_start\tprobe_end\tid\n"
        "chr1\t100\t120\tp1\n",
        encoding="utf-8",
    )
    out_table = tmp_path / "out.tsv"

    ann_processor = probe_module.AnnDataProcessor(ann_table=ann_table)
    config = probe_module.ProcessingConfig(threads=2, variant_cutoff=0, indel_cutoff=0)
    variant_overlap_df = pd.DataFrame({"id": ["p1"], "variant_overlap": [1]})

    probe_module.write_filtered_output_in_chunks(
        ann_data_processor=ann_processor,
        out_table=out_table,
        config=config,
        variant_overlap_df=variant_overlap_df,
        indel_overlap_df=None,
        id_list=None,
    )

    result_df = pd.read_table(out_table)
    assert result_df.empty
    assert list(result_df.columns) == [
        "chrom",
        "probe_start",
        "probe_end",
        "id",
        "indel_overlap",
        "variant_overlap",
    ]


def test_cli_merged_vcf_returns_exit_code_1_on_pipeline_error(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    ann_table = tmp_path / "ann.tsv"
    vcf_path = tmp_path / "test.vcf"
    out_table = tmp_path / "out.tsv"
    ann_table.write_text("chrom\tprobe_start\tprobe_end\tid\n", encoding="utf-8")
    vcf_path.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")

    def fake_pipeline(**kwargs) -> None:
        raise RuntimeError("boom")

    monkeypatch.setattr(probe_module, "run_merged_vcf_pipeline", fake_pipeline)

    result = runner.invoke(
        probe_module.app,
        ["merged-vcf", str(ann_table), str(vcf_path), str(out_table)],
    )

    assert result.exit_code == 1
    assert "处理失败: boom" in result.output
