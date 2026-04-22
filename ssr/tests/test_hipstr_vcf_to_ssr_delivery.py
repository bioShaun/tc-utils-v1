#!/usr/bin/env python3
"""Tests for ssr/hipstr_vcf_to_ssr_delivery.py."""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest
import typer

sys.path.insert(0, str(Path(__file__).parent.parent))

from hipstr_vcf_to_ssr_delivery import (
    BASE_COLUMNS,
    PRIMER_COLUMNS,
    STAT_COLUMNS,
    FilterCriteria,
    ParsedRecord,
    RunSummary,
    annotate_primers,
    annotate_specificity,
    calculate_pic,
    chunked,
    clean_final_dataframe,
    complete_gt,
    configure_logging,
    create_dataframe,
    derive_motif,
    design_primer_for_row,
    ensure_blast_database,
    format_output_dataframe,
    format_sample_genotype,
    get_final_columns,
    import_primer3_module,
    normalize_fraction,
    parse_record,
    parse_vcf_rows,
    passes_initial_filters,
    render_report,
    run_blast_batch,
    run_pipeline,
    split_gb_values,
    to_repeat_values,
    update_summary_from_output,
    version_callback,
    write_outputs,
)


class MockSampleCall:
    def __init__(self, gt: Any, gb: Any) -> None:
        self._gt = gt
        self._gb = gb

    def get(self, key: str) -> Any:
        if key == "GT":
            return self._gt
        if key == "GB":
            return self._gb
        return None


class MockInfo:
    def __init__(self, data: dict[str, Any]) -> None:
        self._data = data

    def get(self, key: str, default: Any = None) -> Any:
        return self._data.get(key, default)


class MockRecord:
    def __init__(
        self,
        contig: str,
        pos: int,
        ref: str,
        stop: int,
        info: dict[str, Any],
        samples: dict[str, MockSampleCall],
    ) -> None:
        self.contig = contig
        self.pos = pos
        self.ref = ref
        self.stop = stop
        self.info = MockInfo(info)
        self.samples = samples


def make_record(
    contig: str = "chr1",
    pos: int = 100,
    ref: str = "ATATAT",
    stop: int = 105,
    start_info: int | None = None,
    end_info: int | None = None,
    period: int = 2,
    samples: dict[str, MockSampleCall] | None = None,
) -> MockRecord:
    info: dict[str, Any] = {"PERIOD": period}
    if start_info is not None:
        info["START"] = start_info
    if end_info is not None:
        info["END"] = end_info
    return MockRecord(
        contig=contig,
        pos=pos,
        ref=ref,
        stop=stop,
        info=info,
        samples=samples or {},
    )


class TestDeriveMotif:
    def test_normal(self) -> None:
        assert derive_motif("ATATAT", pos=100, start_info=100, period=2) == "AT"

    def test_start_before_pos(self) -> None:
        assert derive_motif("ATATAT", pos=100, start_info=98, period=2) == "AT"

    def test_empty_ref(self) -> None:
        assert derive_motif("", pos=100, start_info=100, period=2) == ""

    def test_zero_period(self) -> None:
        assert derive_motif("ATATAT", pos=100, start_info=100, period=0) == ""

    def test_fallback_to_ref_prefix(self) -> None:
        assert derive_motif("ATATAT", pos=100, start_info=200, period=2) == "AT"


class TestSplitGbValues:
    def test_none(self) -> None:
        assert split_gb_values(None) == []

    def test_pipe_separated(self) -> None:
        assert split_gb_values("0|2") == ["0", "2"]

    def test_slash_separated(self) -> None:
        assert split_gb_values("0/2") == ["0", "2"]

    def test_empty_parts_ignored(self) -> None:
        assert split_gb_values("0/") == ["0"]


class TestCompleteGt:
    def test_none(self) -> None:
        assert complete_gt(None) is None

    def test_valid(self) -> None:
        assert complete_gt((0, 1)) == (0, 1)

    def test_with_none_allele(self) -> None:
        assert complete_gt((0, None)) is None

    def test_non_int_allele(self) -> None:
        assert complete_gt((0, ".")) is None


class TestToRepeatValues:
    def test_none_gt(self) -> None:
        call = MockSampleCall(gt=None, gb="0|0")
        assert to_repeat_values(call, ref_repeat=5, period=2) is None

    def test_valid(self) -> None:
        call = MockSampleCall(gt=(0, 1), gb="0|2")
        assert to_repeat_values(call, ref_repeat=5, period=2) == [5, 6]

    def test_missing_gb(self) -> None:
        call = MockSampleCall(gt=(0, 1), gb=None)
        assert to_repeat_values(call, ref_repeat=5, period=2) is None

    def test_dot_in_gb(self) -> None:
        call = MockSampleCall(gt=(0, 1), gb=".|2")
        assert to_repeat_values(call, ref_repeat=5, period=2) is None

    def test_short_gb(self) -> None:
        call = MockSampleCall(gt=(0, 1), gb="0")
        assert to_repeat_values(call, ref_repeat=5, period=2) is None


class TestFormatSampleGenotype:
    def test_empty(self) -> None:
        assert format_sample_genotype([]) == ""

    def test_homozygous(self) -> None:
        assert format_sample_genotype([5, 5]) == "5"

    def test_heterozygous(self) -> None:
        assert format_sample_genotype([5, 7]) == "5/7"


class TestCalculatePic:
    def test_empty(self) -> None:
        assert calculate_pic([]) is None

    def test_single_allele(self) -> None:
        assert calculate_pic([0, 0, 0]) == pytest.approx(0.0, abs=1e-4)

    def test_two_alleles_equal_freq(self) -> None:
        assert calculate_pic([0, 1]) == pytest.approx(0.375, abs=1e-4)

    def test_three_alleles(self) -> None:
        allele_indices = [0, 1, 2]
        freqs = [1 / 3, 1 / 3, 1 / 3]
        sum_pi_sq = sum(f * f for f in freqs)
        second_term = 0.0
        for i, fi in enumerate(freqs):
            for fj in freqs[i + 1 :]:
                second_term += 2 * (fi * fi) * (fj * fj)
        expected = round(1 - sum_pi_sq - second_term, 4)
        assert calculate_pic(allele_indices) == pytest.approx(expected, abs=1e-4)


class TestNormalizeFraction:
    def test_fraction(self) -> None:
        assert normalize_fraction(0.8, "x") == 0.8

    def test_percentage(self) -> None:
        assert normalize_fraction(80, "x") == 0.8

    def test_negative(self) -> None:
        with pytest.raises(typer.BadParameter):
            normalize_fraction(-1, "x")

    def test_too_large(self) -> None:
        with pytest.raises(typer.BadParameter):
            normalize_fraction(101, "x")


class TestChunked:
    def test_empty(self) -> None:
        assert chunked([], size=2) == []

    def test_exact(self) -> None:
        assert chunked([("a", "b", "c")] * 4, size=2) == [
            [("a", "b", "c"), ("a", "b", "c")],
            [("a", "b", "c"), ("a", "b", "c")],
        ]

    def test_remainder(self) -> None:
        assert chunked([("a", "b", "c")] * 3, size=2) == [
            [("a", "b", "c"), ("a", "b", "c")],
            [("a", "b", "c")],
        ]


class TestGetFinalColumns:
    def test_structure(self) -> None:
        samples = ["S1", "S2"]
        cols = get_final_columns(samples)
        assert cols[: len(BASE_COLUMNS)] == BASE_COLUMNS
        assert cols[len(BASE_COLUMNS) : len(BASE_COLUMNS) + len(samples)] == samples
        tail = cols[len(BASE_COLUMNS) + len(samples) :]
        assert tail == STAT_COLUMNS + PRIMER_COLUMNS


class TestCreateDataFrame:
    def test_empty(self) -> None:
        df = create_dataframe([], sample_names=["S1"])
        assert list(df.columns) == get_final_columns(["S1"])
        assert df.empty

    def test_with_rows(self) -> None:
        rows = [{"SSR_ID": "x", "seq_id": "chr1", "start": 1, "end": 6}]
        df = create_dataframe(rows, sample_names=[])
        assert not df.empty
        for col in PRIMER_COLUMNS:
            assert col in df.columns


class TestFormatOutputDataFrame:
    def test_round_and_int(self) -> None:
        rows = [
            {
                "SSR_ID": "x",
                "seq_id": "chr1",
                "start": 1.0,
                "end": 6.0,
                "unit_size": 2.0,
                "ref_repeat": 3.0,
                "maf_like": 0.123456,
                "het_rate": 0.123456,
                "Left_Tm": 60.123456,
                "Right_Tm": 60.123456,
                "missing_count": 0.0,
                "allele_num": 2.0,
                "repeat_min": 3.0,
                "repeat_max": 5.0,
                "repeat_diff": 2.0,
                "bp_diff": 4.0,
                "sample_count": 2.0,
                "Product_Size": 200.0,
            }
        ]
        df = pd.DataFrame(rows)
        out = format_output_dataframe(df, sample_names=[])
        assert out["maf_like"].iloc[0] == 0.123
        assert out["het_rate"].iloc[0] == 0.123
        assert out["start"].iloc[0] == 1
        assert out["Product_Size"].iloc[0] == 200
        assert out.isna().sum().sum() == 0


class TestCleanFinalDataFrame:
    def test_noop_when_disabled(self) -> None:
        df = pd.DataFrame({"repeat_min": [-1, 0, 1]})
        summary = RunSummary(sample_names=[])
        result = clean_final_dataframe(df, summary, clean_final=False, primers_attempted=False)
        assert len(result) == 3

    def test_remove_negative_repeat_min(self) -> None:
        df = pd.DataFrame({"repeat_min": [-1, 0, 1]})
        summary = RunSummary(sample_names=[])
        result = clean_final_dataframe(df, summary, clean_final=True, primers_attempted=False)
        assert len(result) == 2
        assert summary.negative_repeat_removed == 1

    def test_remove_missing_primers(self) -> None:
        df = pd.DataFrame(
            {
                "repeat_min": [1, 1, 1],
                "Left_Primer": ["", "A", "A"],
                "Right_Primer": ["", "", "B"],
            }
        )
        summary = RunSummary(sample_names=[])
        result = clean_final_dataframe(df, summary, clean_final=True, primers_attempted=True)
        assert len(result) == 1
        assert summary.missing_primer_removed == 2


class TestUpdateSummaryFromOutput:
    def test_basic(self) -> None:
        df = pd.DataFrame(
            {
                "unit_size": [2, 3, 2],
                "Left_Primer": ["A", "", "C"],
                "Specificity": ["Specific", "Non-specific", ""],
            }
        )
        summary = RunSummary(sample_names=[])
        update_summary_from_output(df, summary)
        assert summary.final_rows == 3
        assert summary.period_counts == Counter({2: 2, 3: 1})
        assert summary.primer_rows_with_sequences == 2
        assert summary.specificity_counts == Counter({"Specific": 1, "Non-specific": 1})


class TestWriteOutputs:
    def test_creates_files(self, tmp_path: Path) -> None:
        df = pd.DataFrame({"a": [1]})
        tsv = tmp_path / "out.tsv"
        xlsx = tmp_path / "out.xlsx"
        write_outputs(df, tsv, xlsx)
        assert tsv.exists()
        assert xlsx.exists()


class TestRenderReport:
    def test_basic(self, tmp_path: Path) -> None:
        report = tmp_path / "report.txt"
        criteria = FilterCriteria(
            alt_alleles_min=4,
            alt_alleles_max=7,
            pic_min=0.25,
            call_rate_min=0.8,
            skip_initial_filters=False,
        )
        summary = RunSummary(
            sample_names=["S1"],
            total_input_loci=10,
            passed_initial_filters=5,
            final_rows=3,
            period_counts=Counter({2: 2, 3: 1}),
        )
        render_report(
            vcf_path=Path("in.vcf"),
            report_path=report,
            criteria=criteria,
            summary=summary,
            genome=None,
            out_tsv=Path("out.tsv"),
            out_xlsx=Path("out.xlsx"),
        )
        text = report.read_text()
        assert "SSR Filtering Report" in text
        assert "Period 2: 2 loci" in text


class TestVersionCallback:
    def test_raises_exit(self) -> None:
        with pytest.raises(typer.Exit):
            version_callback(True)


class TestConfigureLogging:
    def test_runs(self) -> None:
        configure_logging(verbose=True)
        configure_logging(verbose=False)


class TestParseRecord:
    def test_simple(self) -> None:
        record = make_record(
            pos=100,
            ref="ATATAT",
            stop=105,
            start_info=100,
            end_info=105,
            period=2,
            samples={
                "S1": MockSampleCall(gt=(0, 1), gb="0|2"),
                "S2": MockSampleCall(gt=(0, 0), gb="0|0"),
            },
        )
        parsed = parse_record(record, sample_names=["S1", "S2"])
        assert parsed.row["SSR_ID"] == "SSR_chr1_100"
        assert parsed.row["seq_id"] == "chr1"
        assert parsed.row["start"] == 100
        assert parsed.row["end"] == 105
        assert parsed.row["motif"] == "AT"
        assert parsed.row["unit_size"] == 2
        assert parsed.row["ref_repeat"] == 3
        assert parsed.row["sample_count"] == 2
        assert parsed.row["S1"] == "3/4"
        assert parsed.row["S2"] == "3"
        assert parsed.observed_allele_count == 2
        assert parsed.pic is not None
        assert parsed.call_rate == 1.0

    def test_missing_gt(self) -> None:
        record = make_record(
            samples={
                "S1": MockSampleCall(gt=None, gb="0|0"),
            },
        )
        parsed = parse_record(record, sample_names=["S1"])
        assert parsed.row["S1"] == ""
        assert parsed.row["sample_count"] == 0
        assert parsed.call_rate == 0.0

    def test_invalid_period_raises(self) -> None:
        record = make_record(period=0)
        with pytest.raises(ValueError, match="Invalid PERIOD"):
            parse_record(record, sample_names=[])


class TestPassesInitialFilters:
    def test_skip(self) -> None:
        criteria = FilterCriteria(2, 5, 0.25, 0.8, skip_initial_filters=True)
        summary = RunSummary(sample_names=[])
        parsed = ParsedRecord({}, 0, 0.0, 0.0, 2)
        assert passes_initial_filters(parsed, criteria, summary) is True
        assert summary.passed_initial_filters == 1

    def test_fail_alt(self) -> None:
        criteria = FilterCriteria(2, 5, 0.25, 0.8, skip_initial_filters=False)
        summary = RunSummary(sample_names=[])
        parsed = ParsedRecord({}, observed_allele_count=1, pic=0.5, call_rate=1.0, period=2)
        assert passes_initial_filters(parsed, criteria, summary) is False
        assert summary.failed_alt_filter == 1

    def test_fail_pic(self) -> None:
        criteria = FilterCriteria(2, 5, 0.25, 0.8, skip_initial_filters=False)
        summary = RunSummary(sample_names=[])
        parsed = ParsedRecord({}, observed_allele_count=3, pic=0.1, call_rate=1.0, period=2)
        assert passes_initial_filters(parsed, criteria, summary) is False
        assert summary.failed_pic_filter == 1

    def test_fail_call_rate(self) -> None:
        criteria = FilterCriteria(2, 5, 0.25, 0.8, skip_initial_filters=False)
        summary = RunSummary(sample_names=[])
        parsed = ParsedRecord({}, observed_allele_count=3, pic=0.5, call_rate=0.5, period=2)
        assert passes_initial_filters(parsed, criteria, summary) is False
        assert summary.failed_call_rate_filter == 1


class TestParseVcfRows:
    def test_basic(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "test.vcf"
        vcf_path.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=1000>\n"
            "##INFO=<ID=PERIOD,Number=1,Type=Integer,Description=\"Period\">\n"
            "##INFO=<ID=START,Number=1,Type=Integer,Description=\"Start\">\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            "##FORMAT=<ID=GB,Number=1,Type=String,Description=\"GB\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr1\t100\t.\tATATAT\tA\t.\t.\tPERIOD=2;START=100;END=105\tGT:GB\t0/1:0|2\n"
        )
        criteria = FilterCriteria(1, 10, 0.0, 0.0, skip_initial_filters=True)
        sample_names, rows, summary = parse_vcf_rows(vcf_path, criteria)
        assert sample_names == ["S1"]
        assert len(rows) == 1
        assert summary.total_input_loci == 1

    def test_filters_monomorphic(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "test.vcf"
        vcf_path.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=1000>\n"
            "##INFO=<ID=PERIOD,Number=1,Type=Integer,Description=\"Period\">\n"
            "##INFO=<ID=START,Number=1,Type=Integer,Description=\"Start\">\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            "##FORMAT=<ID=GB,Number=1,Type=String,Description=\"GB\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr1\t100\t.\tATATAT\t.\t.\t.\tPERIOD=2;START=100;END=105\tGT:GB\t0/0:0|0\n"
        )
        criteria = FilterCriteria(1, 10, 0.0, 0.0, skip_initial_filters=True)
        sample_names, rows, summary = parse_vcf_rows(vcf_path, criteria)
        assert len(rows) == 0
        assert summary.monomorphic_removed == 1


class MockFasta:
    def __init__(self, seqs: dict[str, str]) -> None:
        self._seqs = seqs

    def get_reference_length(self, chrom: str) -> int:
        return len(self._seqs[chrom])

    def fetch(self, chrom: str, start: int, end: int) -> str:
        return self._seqs[chrom][start:end]


class TestDesignPrimerForRow:
    def test_no_primer_when_coords_bad(self) -> None:
        row = pd.Series({"SSR_ID": "x", "seq_id": "chr1", "start": 0, "end": 5})
        fasta = MockFasta({"chr1": "A" * 1000})
        result = design_primer_for_row(row, fasta=fasta, primer3=MagicMock())
        assert result["Left_Primer"] == ""

    def test_no_primer_when_contig_missing(self) -> None:
        row = pd.Series({"SSR_ID": "x", "seq_id": "chr2", "start": 10, "end": 15})
        fasta = MockFasta({"chr1": "A" * 1000})
        result = design_primer_for_row(row, fasta=fasta, primer3=MagicMock())
        assert result["Left_Primer"] == ""

    def test_success(self) -> None:
        row = pd.Series({"SSR_ID": "x", "seq_id": "chr1", "start": 401, "end": 406})
        seq = "A" * 1000
        fasta = MockFasta({"chr1": seq})
        primer3 = MagicMock()
        primer3.bindings.designPrimers.return_value = {
            "PRIMER_LEFT_0_SEQUENCE": "ACGT",
            "PRIMER_RIGHT_0_SEQUENCE": "TGCA",
            "PRIMER_LEFT_0_TM": 60.123,
            "PRIMER_RIGHT_0_TM": 59.987,
            "PRIMER_PAIR_0_PRODUCT_SIZE": 250,
        }
        result = design_primer_for_row(row, fasta=fasta, primer3=primer3)
        assert result["Left_Primer"] == "ACGT"
        assert result["Right_Primer"] == "TGCA"
        assert result["Left_Tm"] == 60.123
        assert result["Right_Tm"] == 59.987
        assert result["Product_Size"] == 250

    def test_fallback_design_primers(self) -> None:
        row = pd.Series({"SSR_ID": "x", "seq_id": "chr1", "start": 401, "end": 406})
        seq = "A" * 1000
        fasta = MockFasta({"chr1": seq})
        primer3 = MagicMock()
        del primer3.bindings.designPrimers
        primer3.bindings.design_primers.return_value = {
            "PRIMER_LEFT_0_SEQUENCE": "ACGT",
            "PRIMER_RIGHT_0_SEQUENCE": "TGCA",
            "PRIMER_LEFT_0_TM": 60.0,
            "PRIMER_RIGHT_0_TM": 60.0,
            "PRIMER_PAIR_0_PRODUCT_SIZE": 200,
        }
        result = design_primer_for_row(row, fasta=fasta, primer3=primer3)
        assert result["Left_Primer"] == "ACGT"

    def test_no_result_when_primer3_missing_keys(self) -> None:
        row = pd.Series({"SSR_ID": "x", "seq_id": "chr1", "start": 401, "end": 406})
        seq = "A" * 1000
        fasta = MockFasta({"chr1": seq})
        primer3 = MagicMock()
        primer3.bindings.designPrimers.return_value = {}
        result = design_primer_for_row(row, fasta=fasta, primer3=primer3)
        assert result["Left_Primer"] == ""

    def test_no_result_when_template_too_short(self) -> None:
        row = pd.Series({"SSR_ID": "x", "seq_id": "chr1", "start": 2, "end": 3})
        seq = "A" * 10
        fasta = MockFasta({"chr1": seq})
        result = design_primer_for_row(row, fasta=fasta, primer3=MagicMock())
        assert result["Left_Primer"] == ""


class TestAnnotatePrimers:
    def test_empty_df(self) -> None:
        df = pd.DataFrame()
        result = annotate_primers(df, genome=Path("/dev/null"))
        assert result.empty

    @patch("hipstr_vcf_to_ssr_delivery.import_primer3_module")
    @patch("hipstr_vcf_to_ssr_delivery.pysam.FastaFile")
    def test_annotates(self, mock_fasta_class, mock_import_primer3) -> None:
        mock_fasta = MagicMock()
        mock_fasta.__enter__ = MagicMock(return_value=mock_fasta)
        mock_fasta.__exit__ = MagicMock(return_value=False)
        mock_fasta.get_reference_length.return_value = 1000
        mock_fasta.fetch.return_value = "A" * 400
        mock_fasta_class.return_value = mock_fasta

        primer3 = MagicMock()
        primer3.bindings.designPrimers.return_value = {
            "PRIMER_LEFT_0_SEQUENCE": "ACGT",
            "PRIMER_RIGHT_0_SEQUENCE": "TGCA",
            "PRIMER_LEFT_0_TM": 60.0,
            "PRIMER_RIGHT_0_TM": 60.0,
            "PRIMER_PAIR_0_PRODUCT_SIZE": 200,
        }
        mock_import_primer3.return_value = primer3

        df = pd.DataFrame(
            {
                "SSR_ID": ["SSR_chr1_100"],
                "seq_id": ["chr1"],
                "start": [401],
                "end": [406],
            }
        )
        result = annotate_primers(df, genome=Path("/fake.fa"))
        assert result["Left_Primer"].iloc[0] == "ACGT"
        assert result["Right_Primer"].iloc[0] == "TGCA"


class TestImportPrimer3Module:
    def test_success(self) -> None:
        fake_module = MagicMock()
        with patch("builtins.__import__", return_value=fake_module):
            mod = import_primer3_module()
            assert mod is fake_module

    def test_failure(self) -> None:
        with patch("builtins.__import__", side_effect=ImportError("no module")):
            with pytest.raises(typer.BadParameter):
                import_primer3_module()


class TestEnsureBlastDatabase:
    def test_reuses_existing(self, tmp_path: Path) -> None:
        genome = tmp_path / "ref.fa"
        genome.write_text(">chr1\nACGT\n")
        db_prefix = tmp_path / "blastdb" / "ref"
        db_prefix.parent.mkdir(parents=True)
        for ext in [".nhr", ".nin", ".nsq"]:
            db_prefix.with_suffix(ext).touch()
        result = ensure_blast_database(genome, tmp_path, "makeblastdb")
        assert result == db_prefix

    @patch("hipstr_vcf_to_ssr_delivery.subprocess.run")
    def test_builds_new(self, mock_run, tmp_path: Path) -> None:
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
        genome = tmp_path / "ref.fa"
        genome.write_text(">chr1\nACGT\n")
        result = ensure_blast_database(genome, tmp_path, "makeblastdb")
        assert result.name == "ref"
        mock_run.assert_called_once()

    @patch("hipstr_vcf_to_ssr_delivery.subprocess.run")
    def test_failure(self, mock_run, tmp_path: Path) -> None:
        mock_run.return_value = MagicMock(returncode=1, stdout="err", stderr="fail")
        genome = tmp_path / "ref.fa"
        genome.write_text(">chr1\nACGT\n")
        with pytest.raises(RuntimeError, match="makeblastdb failed"):
            ensure_blast_database(genome, tmp_path, "makeblastdb")

    def test_binary_not_found(self, tmp_path: Path) -> None:
        genome = tmp_path / "ref.fa"
        genome.write_text(">chr1\nACGT\n")
        with pytest.raises(FileNotFoundError, match="makeblastdb not found"):
            ensure_blast_database(genome, tmp_path, "/nonexistent/makeblastdb")


class TestRunBlastBatch:
    @patch("hipstr_vcf_to_ssr_delivery.shutil.which")
    @patch("hipstr_vcf_to_ssr_delivery.subprocess.run")
    def test_specific(self, mock_run, mock_which, tmp_path: Path) -> None:
        mock_which.return_value = "/usr/bin/blastn"
        stdout = (
            "SSR_1_L\tchr1\t100.0\t8\t0\t0\t1\t8\t1\t8\t1e-10\t40\n"
            "SSR_1_R\tchr1\t100.0\t8\t0\t0\t1\t8\t1\t8\t1e-10\t40\n"
        )
        mock_run.return_value = MagicMock(returncode=0, stdout=stdout, stderr="")
        rows = [("SSR_1", "ACGTACGT", "TGCATGCA")]
        result = run_blast_batch(rows, blastn_bin="blastn", blast_db=tmp_path / "db")
        assert result["SSR_1"] == "Specific"

    @patch("hipstr_vcf_to_ssr_delivery.shutil.which")
    @patch("hipstr_vcf_to_ssr_delivery.subprocess.run")
    def test_non_specific(self, mock_run, mock_which, tmp_path: Path) -> None:
        mock_which.return_value = "/usr/bin/blastn"
        stdout = (
            "SSR_1_L\tchr1\t100.0\t8\t0\t0\t1\t8\t1\t8\t1e-10\t40\n"
            "SSR_1_R\tchr1\t100.0\t8\t0\t0\t1\t8\t1\t8\t1e-10\t40\n"
            "SSR_1_R\tchr2\t100.0\t8\t0\t0\t1\t8\t1\t8\t1e-10\t40\n"
        )
        mock_run.return_value = MagicMock(returncode=0, stdout=stdout, stderr="")
        rows = [("SSR_1", "ACGTACGT", "TGCATGCA")]
        result = run_blast_batch(rows, blastn_bin="blastn", blast_db=tmp_path / "db")
        assert "Non-specific" in result["SSR_1"]
        assert "L:1" in result["SSR_1"]
        assert "R:2" in result["SSR_1"]

    @patch("hipstr_vcf_to_ssr_delivery.shutil.which")
    @patch("hipstr_vcf_to_ssr_delivery.subprocess.run")
    def test_blastn_failure(self, mock_run, mock_which, tmp_path: Path) -> None:
        mock_which.return_value = "/usr/bin/blastn"
        mock_run.return_value = MagicMock(returncode=1, stdout="", stderr="error")
        rows = [("SSR_1", "ACGT", "TGCA")]
        with pytest.raises(RuntimeError, match="blastn failed"):
            run_blast_batch(rows, blastn_bin="blastn", blast_db=tmp_path / "db")

    def test_blastn_not_found(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="blastn not found"):
            run_blast_batch([], blastn_bin="/nonexistent/blastn", blast_db=tmp_path / "db")


class TestAnnotateSpecificity:
    def test_empty_df(self) -> None:
        df = pd.DataFrame()
        result = annotate_specificity(df, blast_db=Path("db"), blastn_bin="blastn", batch_size=10)
        assert result.empty

    def test_no_primers(self) -> None:
        df = pd.DataFrame(
            {
                "SSR_ID": ["SSR_1"],
                "Left_Primer": [""],
                "Right_Primer": [""],
            }
        )
        result = annotate_specificity(df, blast_db=Path("db"), blastn_bin="blastn", batch_size=10)
        assert "Specificity" not in result.columns

    @patch("hipstr_vcf_to_ssr_delivery.run_blast_batch")
    @patch("hipstr_vcf_to_ssr_delivery.shutil.which")
    def test_with_primers(self, mock_which, mock_batch) -> None:
        mock_which.return_value = "/usr/bin/blastn"
        mock_batch.return_value = {"SSR_1": "Specific"}
        df = pd.DataFrame(
            {
                "SSR_ID": ["SSR_1"],
                "Left_Primer": ["ACGT"],
                "Right_Primer": ["TGCA"],
            }
        )
        result = annotate_specificity(df, blast_db=Path("db"), blastn_bin="blastn", batch_size=10)
        assert result["Specificity"].iloc[0] == "Specific"


class TestRunPipeline:
    def test_minimal(self, tmp_path: Path) -> None:
        vcf = tmp_path / "in.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
            "chr1\t100\t.\tATATAT\t.\t.\t.\tPERIOD=2;START=100;END=105\tGT:GB\t0/1:0|2\t0/0:0|0\n"
        )
        out_xlsx = tmp_path / "out.xlsx"
        run_pipeline(
            vcf=vcf,
            out_xlsx=out_xlsx,
            genome=None,
            blast_db=None,
            workdir=None,
            out_tsv=None,
            report_txt=None,
            alleles_min=1,
            alleles_max=10,
            pic_min=0.0,
            call_rate_min=0.0,
            skip_initial_filters=True,
            clean_final=True,
            check_specificity=False,
            batch_size=10,
            blastn="blastn",
            makeblastdb="makeblastdb",
        )
        assert out_xlsx.exists()
        tsv = out_xlsx.with_suffix(".tsv")
        assert tsv.exists()
        report = out_xlsx.with_name("filtering_report.txt")
        assert report.exists()

    @patch("hipstr_vcf_to_ssr_delivery.pysam.FastaFile")
    @patch("hipstr_vcf_to_ssr_delivery.import_primer3_module")
    def test_with_genome(self, mock_import_primer3, mock_fasta_class, tmp_path: Path) -> None:
        vcf = tmp_path / "in.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=1000>\n"
            "##INFO=<ID=PERIOD,Number=1,Type=Integer,Description=\"Period\">\n"
            "##INFO=<ID=START,Number=1,Type=Integer,Description=\"Start\">\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            "##FORMAT=<ID=GB,Number=1,Type=String,Description=\"GB\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr1\t100\t.\tATATAT\tA\t.\t.\tPERIOD=2;START=100;END=105\tGT:GB\t0/1:0|2\n"
        )
        genome = tmp_path / "genome.fa"
        genome.write_text(">chr1\n" + "A" * 1000 + "\n")

        mock_fasta = MagicMock()
        mock_fasta.__enter__ = MagicMock(return_value=mock_fasta)
        mock_fasta.__exit__ = MagicMock(return_value=False)
        mock_fasta.get_reference_length.return_value = 1000
        mock_fasta.fetch.return_value = "A" * 400
        mock_fasta_class.return_value = mock_fasta

        primer3 = MagicMock()
        primer3.bindings.designPrimers.return_value = {
            "PRIMER_LEFT_0_SEQUENCE": "ACGT",
            "PRIMER_RIGHT_0_SEQUENCE": "TGCA",
            "PRIMER_LEFT_0_TM": 60.0,
            "PRIMER_RIGHT_0_TM": 60.0,
            "PRIMER_PAIR_0_PRODUCT_SIZE": 200,
        }
        mock_import_primer3.return_value = primer3

        out_xlsx = tmp_path / "out.xlsx"
        run_pipeline(
            vcf=vcf,
            out_xlsx=out_xlsx,
            genome=genome,
            blast_db=None,
            workdir=None,
            out_tsv=None,
            report_txt=None,
            alleles_min=1,
            alleles_max=10,
            pic_min=0.0,
            call_rate_min=0.0,
            skip_initial_filters=True,
            clean_final=True,
            check_specificity=False,
            batch_size=10,
            blastn="blastn",
            makeblastdb="makeblastdb",
        )
        assert out_xlsx.exists()
        df = pd.read_excel(out_xlsx)
        assert df["Left_Primer"].iloc[0] == "ACGT"


class TestMainCli:
    @patch("hipstr_vcf_to_ssr_delivery.run_pipeline")
    def test_alleles_min_max_validation(self, mock_run, tmp_path: Path) -> None:
        vcf = tmp_path / "in.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        out_xlsx = tmp_path / "out.xlsx"
        with pytest.raises(typer.BadParameter, match="cannot be greater"):
            from hipstr_vcf_to_ssr_delivery import main
            main(
                vcf=vcf,
                out_xlsx=out_xlsx,
                genome=None,
                blast_db=None,
                workdir=None,
                out_tsv=None,
                report_txt=None,
                alleles_min=10,
                alleles_max=5,
                pic_min=0.25,
                call_rate_min=0.8,
                skip_initial_filters=False,
                clean_final=True,
                check_specificity=True,
                batch_size=100,
                blastn="blastn",
                makeblastdb="makeblastdb",
                verbose=False,
                version=None,
            )

    @patch("hipstr_vcf_to_ssr_delivery.run_pipeline")
    def test_batch_size_validation(self, mock_run, tmp_path: Path) -> None:
        vcf = tmp_path / "in.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        out_xlsx = tmp_path / "out.xlsx"
        with pytest.raises(typer.BadParameter, match="positive integer"):
            from hipstr_vcf_to_ssr_delivery import main
            main(
                vcf=vcf,
                out_xlsx=out_xlsx,
                genome=None,
                blast_db=None,
                workdir=None,
                out_tsv=None,
                report_txt=None,
                alleles_min=1,
                alleles_max=5,
                pic_min=0.25,
                call_rate_min=0.8,
                skip_initial_filters=False,
                clean_final=True,
                check_specificity=True,
                batch_size=0,
                blastn="blastn",
                makeblastdb="makeblastdb",
                verbose=False,
                version=None,
            )

    @patch("hipstr_vcf_to_ssr_delivery.run_pipeline")
    def test_genome_not_found(self, mock_run, tmp_path: Path) -> None:
        vcf = tmp_path / "in.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        out_xlsx = tmp_path / "out.xlsx"
        with pytest.raises(typer.BadParameter, match="Genome FASTA not found"):
            from hipstr_vcf_to_ssr_delivery import main
            main(
                vcf=vcf,
                out_xlsx=out_xlsx,
                genome=Path("/nonexistent.fa"),
                blast_db=None,
                workdir=None,
                out_tsv=None,
                report_txt=None,
                alleles_min=1,
                alleles_max=5,
                pic_min=0.25,
                call_rate_min=0.8,
                skip_initial_filters=False,
                clean_final=True,
                check_specificity=True,
                batch_size=100,
                blastn="blastn",
                makeblastdb="makeblastdb",
                verbose=False,
                version=None,
            )
