"""Tests for primer/kasp_mapper.py.

Covers:
- parse_ssr (3-column SSR TSV parsing)
- detect_format (format1 / format2 / ssr discrimination)
- load_id_chr_map (ID -> target chromosomes hints)
- _find_kasp_snp_query_positions (3'-based SNP finder)
- analyze_format1_loci (subject SNP via _query_to_subject_pos on both FAM & HEX)
- analyze_format2_loci (with target_chroms)
- analyze_ssr_loci (strand / distance / coverage / amplicon filters)
- create_blast_query_file (SSR branch writes _F/_R entries)
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pytest

PRIMER_DIR = Path(__file__).resolve().parent.parent
if str(PRIMER_DIR) not in sys.path:
    sys.path.insert(0, str(PRIMER_DIR))

kasp_mapper = importlib.import_module("kasp_mapper")

BlastResult = kasp_mapper.BlastResult
KaspPrimer = kasp_mapper.KaspPrimer


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def make_hit(
    query_id: str,
    subject_id: str,
    query_start: int,
    query_end: int,
    subject_start: int,
    subject_end: int,
    *,
    identity: float = 99.0,
    alignment_length: int | None = None,
    mismatches: int = 0,
    gap_opens: int = 0,
    evalue: float = 1e-20,
    bit_score: float = 100.0,
) -> BlastResult:
    if alignment_length is None:
        alignment_length = abs(query_end - query_start) + 1
    return BlastResult(
        query_id=query_id,
        subject_id=subject_id,
        identity=identity,
        alignment_length=alignment_length,
        mismatches=mismatches,
        gap_opens=gap_opens,
        query_start=query_start,
        query_end=query_end,
        subject_start=subject_start,
        subject_end=subject_end,
        evalue=evalue,
        bit_score=bit_score,
    )


# ---------------------------------------------------------------------------
# parse_ssr
# ---------------------------------------------------------------------------
class TestParseSsr:
    def test_parses_three_column_tsv(self, tmp_path: Path) -> None:
        tsv = tmp_path / "ssr.tsv"
        tsv.write_text(
            "SSR001\tACGTACGTACGT\tTTTTCCCCAAAA\n"
            "SSR002\tGGGGAAAACCCC\tAAAATTTTGGGG\n",
            encoding="utf-8",
        )
        primers = kasp_mapper.parse_ssr(tsv)
        assert len(primers) == 2
        assert primers[0].primer_id == "SSR001"
        assert primers[0].format_type == "ssr"
        assert primers[0].forward_primer == "ACGTACGTACGT"
        assert primers[0].reverse_primer == "TTTTCCCCAAAA"

    def test_skips_header_row(self, tmp_path: Path) -> None:
        tsv = tmp_path / "ssr.tsv"
        tsv.write_text(
            "id\tforward\treverse\n"
            "SSR001\tACGTACGT\tTTTTCCCC\n",
            encoding="utf-8",
        )
        primers = kasp_mapper.parse_ssr(tsv)
        assert [p.primer_id for p in primers] == ["SSR001"]

    def test_skips_non_dna_records(self, tmp_path: Path) -> None:
        tsv = tmp_path / "ssr.tsv"
        tsv.write_text(
            "SSR001\tACGTACGT\tTTTTCCCC\n"
            "SSR002\tACGTXXXX\tTTTTCCCC\n",
            encoding="utf-8",
        )
        primers = kasp_mapper.parse_ssr(tsv)
        assert [p.primer_id for p in primers] == ["SSR001"]

    def test_empty_file_raises_value_error(self, tmp_path: Path) -> None:
        tsv = tmp_path / "empty.tsv"
        tsv.write_text("", encoding="utf-8")
        with pytest.raises(ValueError):
            kasp_mapper.parse_ssr(tsv)


# ---------------------------------------------------------------------------
# detect_format
# ---------------------------------------------------------------------------
class TestDetectFormat:
    def test_detects_ssr_for_pure_atgc_three_columns(self, tmp_path: Path) -> None:
        tsv = tmp_path / "ssr.tsv"
        tsv.write_text(
            "SSR001\tACGTACGTACGT\tTTTTCCCCAAAA\n"
            "SSR002\tGGGGAAAACCCC\tAAAATTTTGGGG\n",
            encoding="utf-8",
        )
        assert kasp_mapper.detect_format(tsv) == "ssr"

    def test_detects_format1_when_dye_tags_present(self, tmp_path: Path) -> None:
        tsv = tmp_path / "f1.tsv"
        tsv.write_text(
            f"KASP1\t{kasp_mapper.FAM_TAG}ACGTACGT\t{kasp_mapper.HEX_TAG}ACGTACGT\tTTTTCCCCAAAA\n",
            encoding="utf-8",
        )
        assert kasp_mapper.detect_format(tsv) == "format1"

    def test_detects_format2_when_snp_marker_present(self, tmp_path: Path) -> None:
        tsv = tmp_path / "f2.tsv"
        tsv.write_text("KASP1\tACGT[A/T]CGTACGTACGT\n", encoding="utf-8")
        assert kasp_mapper.detect_format(tsv) == "format2"


# ---------------------------------------------------------------------------
# create_blast_query_file
# ---------------------------------------------------------------------------
class TestCreateBlastQueryFileSsr:
    def test_writes_forward_and_reverse_records(self, tmp_path: Path) -> None:
        primers = [
            KaspPrimer(
                primer_id="SSR001",
                format_type="ssr",
                forward_primer="ACGTACGT",
                reverse_primer="TTTTCCCC",
            )
        ]
        out = tmp_path / "query.fa"
        kasp_mapper.create_blast_query_file(primers, out)
        content = out.read_text(encoding="utf-8")
        assert ">SSR001_F\nACGTACGT\n" in content
        assert ">SSR001_R\nTTTTCCCC\n" in content


# ---------------------------------------------------------------------------
# load_id_chr_map
# ---------------------------------------------------------------------------
class TestLoadIdChrMap:
    def test_returns_empty_dict_when_path_is_none(self) -> None:
        assert kasp_mapper.load_id_chr_map(None) == {}

    def test_parses_tab_separated_two_columns(self, tmp_path: Path) -> None:
        f = tmp_path / "map.tsv"
        f.write_text(
            "primer_id\tchrom\n"
            "KASP1\tchr1\n"
            "KASP2\tchr2,chr3\n"
            "# comment\n"
            "KASP3\tchr4|chr5\n",
            encoding="utf-8",
        )
        mapping = kasp_mapper.load_id_chr_map(f)
        assert mapping["KASP1"] == {"chr1"}
        assert mapping["KASP2"] == {"chr2", "chr3"}
        assert mapping["KASP3"] == {"chr4", "chr5"}

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            kasp_mapper.load_id_chr_map(tmp_path / "missing.tsv")


# ---------------------------------------------------------------------------
# _find_kasp_snp_query_positions (3'-based SNP finder)
# ---------------------------------------------------------------------------
class TestFindKaspSnpQueryPositions:
    def test_finds_3prime_most_mismatch(self) -> None:
        # Mismatches at positions 3 (A/T, from 5') and 8 (3'-most) — biological
        # KASP SNP is at the 3' end; the 5' mismatch is a tuning artifact.
        fam = "ACATTTGA"   # length 8
        hex_ = "ACTTTTGT"  # differs at index 2 and index 7 (0-based)
        result = kasp_mapper._find_kasp_snp_query_positions(fam, hex_)
        assert result is not None
        fam_qpos, hex_qpos, fam_allele, hex_allele = result
        assert fam_qpos == 8 and hex_qpos == 8
        assert fam_allele == "A" and hex_allele == "T"

    def test_handles_different_lengths_using_3prime_anchor(self) -> None:
        # FAM has extra 5' bases; 3'-end comparison should still work.
        fam = "GGGGACGTACA"  # length 11
        hex_ = "ACGTACG"       # length 7
        # 3'-end comparison of last 7 bases: "CGTACA" vs "CGTACG"
        result = kasp_mapper._find_kasp_snp_query_positions(fam, hex_)
        assert result is not None
        fam_qpos, hex_qpos, fam_allele, hex_allele = result
        # 3'-most mismatch: fam last base 'A', hex last base 'G'
        assert fam_qpos == 11
        assert hex_qpos == 7
        assert fam_allele == "A" and hex_allele == "G"

    def test_returns_none_when_identical(self) -> None:
        assert kasp_mapper._find_kasp_snp_query_positions("ACGTACGT", "ACGTACGT") is None


# ---------------------------------------------------------------------------
# analyze_format1_loci
# ---------------------------------------------------------------------------
class TestAnalyzeFormat1Loci:
    def _make_primer(self) -> KaspPrimer:
        # Length 20; mismatch only at the last position -> biological SNP at 3' end.
        fam = "ACGTACGTACGTACGTACGA"
        hex_ = "ACGTACGTACGTACGTACGT"
        return KaspPrimer(
            primer_id="KASP1",
            format_type="format1",
            fam_primer=fam,
            hex_primer=hex_,
            common_primer="GGGGCCCCAAAATTTTGGGG",
        )

    def test_uses_3prime_snp_and_maps_to_subject(self) -> None:
        primer = self._make_primer()
        # FAM maps fully: query 1..20 -> subject 100..119 (forward)
        fam_hit = make_hit("KASP1_FAM", "chr1", 1, 20, 100, 119)
        # HEX maps fully: query 1..20 -> subject 100..119 (forward)
        hex_hit = make_hit("KASP1_HEX", "chr1", 1, 20, 100, 119)
        # Common maps: query 1..20 -> subject 200..219 (forward, downstream)
        common_hit = make_hit("KASP1_Common", "chr1", 1, 20, 200, 219)

        loci = kasp_mapper.analyze_format1_loci(
            primer, [fam_hit, hex_hit, common_hit], min_coverage=0.9
        )
        assert len(loci) == 1
        # 3'-end SNP: query position 20 -> subject 119.
        assert loci[0].snp_pos == 119
        assert loci[0].chrom == "chr1"
        assert loci[0].ref_allele == "A"
        assert loci[0].alt_allele == "T"

    def test_rejects_when_fam_and_hex_snp_positions_disagree(self) -> None:
        primer = self._make_primer()
        # FAM 3' at subject 119 but HEX 3' at subject 200 (100bp apart) -> reject.
        fam_hit = make_hit("KASP1_FAM", "chr1", 1, 20, 100, 119)
        hex_hit = make_hit("KASP1_HEX", "chr1", 1, 20, 181, 200)
        common_hit = make_hit("KASP1_Common", "chr1", 1, 20, 300, 319)
        loci = kasp_mapper.analyze_format1_loci(
            primer, [fam_hit, hex_hit, common_hit], max_snp_distance=10
        )
        assert loci == []

    def test_picks_3prime_snp_even_with_5prime_tuning_mismatch(self) -> None:
        # FAM vs HEX differ at index 2 (5' tuning mismatch) AND at the 3' end.
        # Old logic (first 5' mismatch) -> subject 102. Correct (3'-end) -> 119.
        fam = "ACATACGTACGTACGTACGA"
        hex_ = "ACTTACGTACGTACGTACGT"
        primer = KaspPrimer(
            primer_id="KASP1",
            format_type="format1",
            fam_primer=fam,
            hex_primer=hex_,
            common_primer="GGGGCCCCAAAATTTTGGGG",
        )
        fam_hit = make_hit("KASP1_FAM", "chr1", 1, 20, 100, 119)
        hex_hit = make_hit("KASP1_HEX", "chr1", 1, 20, 100, 119)
        common_hit = make_hit("KASP1_Common", "chr1", 1, 20, 200, 219)
        loci = kasp_mapper.analyze_format1_loci(
            primer, [fam_hit, hex_hit, common_hit], min_coverage=0.9
        )
        assert len(loci) == 1
        # 3'-end SNP, not the tuning 5' mismatch.
        assert loci[0].snp_pos == 119
        assert loci[0].ref_allele == "A"
        assert loci[0].alt_allele == "T"

    def test_target_chroms_filters_hits(self) -> None:
        primer = self._make_primer()
        fam_hit_a = make_hit("KASP1_FAM", "chr1", 1, 20, 100, 119)
        hex_hit_a = make_hit("KASP1_HEX", "chr1", 1, 20, 100, 119)
        common_hit_a = make_hit("KASP1_Common", "chr1", 1, 20, 200, 219)
        fam_hit_b = make_hit("KASP1_FAM", "chr2", 1, 20, 100, 119)
        hex_hit_b = make_hit("KASP1_HEX", "chr2", 1, 20, 100, 119)
        common_hit_b = make_hit("KASP1_Common", "chr2", 1, 20, 200, 219)

        hits = [fam_hit_a, hex_hit_a, common_hit_a, fam_hit_b, hex_hit_b, common_hit_b]
        loci = kasp_mapper.analyze_format1_loci(
            primer, hits, target_chroms={"chr2"}
        )
        assert {loc.chrom for loc in loci} == {"chr2"}


# ---------------------------------------------------------------------------
# analyze_format2_loci
# ---------------------------------------------------------------------------
class TestAnalyzeFormat2Loci:
    def test_target_chroms_narrows_results(self) -> None:
        primer = KaspPrimer(
            primer_id="KASP2",
            format_type="format2",
            flank_seq="A" * 60,
            snp_pos=30,
            ref_allele="A",
            alt_allele="T",
        )
        hit_a = make_hit("KASP2", "chr1", 1, 60, 100, 159)
        hit_b = make_hit("KASP2", "chr2", 1, 60, 300, 359)
        loci, reason = kasp_mapper.analyze_format2_loci(
            primer, [hit_a, hit_b], target_chroms={"chr2"}
        )
        assert reason == "成功"
        assert {loc.chrom for loc in loci} == {"chr2"}

    def test_target_chroms_no_match_reports_reason(self) -> None:
        primer = KaspPrimer(
            primer_id="KASP2",
            format_type="format2",
            flank_seq="A" * 60,
            snp_pos=30,
            ref_allele="A",
            alt_allele="T",
        )
        hit_a = make_hit("KASP2", "chr1", 1, 60, 100, 159)
        loci, reason = kasp_mapper.analyze_format2_loci(
            primer, [hit_a], target_chroms={"chrX"}
        )
        assert loci == []
        assert "chrX" in reason


# ---------------------------------------------------------------------------
# analyze_ssr_loci
# ---------------------------------------------------------------------------
class TestAnalyzeSsrLoci:
    def _make_primer(self) -> KaspPrimer:
        return KaspPrimer(
            primer_id="SSR1",
            format_type="ssr",
            forward_primer="ACGTACGTACGT",      # len 12
            reverse_primer="TTTTCCCCAAAA",      # len 12
        )

    def test_valid_pair_on_opposite_strands_passes(self) -> None:
        primer = self._make_primer()
        # Forward: query 1..12 -> subject 100..111 (+ strand)
        forward_hit = make_hit("SSR1_F", "chr1", 1, 12, 100, 111)
        # Reverse: query 1..12 -> subject 311..300 (- strand, downstream)
        reverse_hit = make_hit("SSR1_R", "chr1", 1, 12, 311, 300)
        loci, reason = kasp_mapper.analyze_ssr_loci(
            primer, [forward_hit, reverse_hit], max_amplicon_size=500
        )
        assert reason == "成功"
        assert len(loci) == 1
        assert loci[0].forward_strand == "+"
        assert loci[0].reverse_strand == "-"
        assert loci[0].amplicon_size == 311 - 100

    def test_same_strand_rejected(self) -> None:
        primer = self._make_primer()
        forward_hit = make_hit("SSR1_F", "chr1", 1, 12, 100, 111)
        reverse_hit = make_hit("SSR1_R", "chr1", 1, 12, 200, 211)  # also + strand
        loci, reason = kasp_mapper.analyze_ssr_loci(
            primer, [forward_hit, reverse_hit]
        )
        assert loci == []
        assert "配对" in reason or "strand" in reason.lower()

    def test_amplicon_too_large_rejected(self) -> None:
        primer = self._make_primer()
        forward_hit = make_hit("SSR1_F", "chr1", 1, 12, 100, 111)
        reverse_hit = make_hit("SSR1_R", "chr1", 1, 12, 5000, 4989)
        loci, reason = kasp_mapper.analyze_ssr_loci(
            primer, [forward_hit, reverse_hit], max_amplicon_size=1000
        )
        assert loci == []
        assert "扩增子" in reason or "配对" in reason

    def test_coverage_below_threshold_rejected(self) -> None:
        primer = self._make_primer()
        # Only 6/12 bases aligned -> coverage 50%
        forward_hit = make_hit("SSR1_F", "chr1", 1, 6, 100, 105)
        reverse_hit = make_hit("SSR1_R", "chr1", 1, 12, 311, 300)
        loci, reason = kasp_mapper.analyze_ssr_loci(
            primer, [forward_hit, reverse_hit], min_coverage=0.9
        )
        assert loci == []
        assert "Coverage" in reason
