from pathlib import Path

import pandas as pd
import pytest

from panel import realign_blast


def _alignment_row(
    *,
    probe_id: str = "probe1",
    chrom: str = "chr1",
    bitscore: int = 100,
    gap_opens: int = 0,
    btop: str | object = "10",
    sstart: int = 100,
    send: int = 109,
    align_len: int = 10,
) -> dict[str, object]:
    strand = "+" if sstart <= send else "-"
    return {
        "id": probe_id,
        "query_len": 10,
        "query_start": 0,
        "strand": strand,
        "chrom": chrom,
        "hit_start": min(sstart, send) - 1,
        "align_len": align_len,
        "bitscore": bitscore,
        "mismatches": 0,
        "gap_opens": gap_opens,
        "qstart": 1,
        "qend": 10,
        "sstart": sstart,
        "send": send,
        "btop": btop,
    }


def _offsets(offset: int = 5) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "id": ["probe1"],
            "offset_fwd": [offset],
            "offset_rev": [10 - offset - 1],
            "alleles": ["A/G"],
        }
    )


def test_fasta_and_offsets_from_probe_table_expands_single_iupac_marker(tmp_path: Path) -> None:
    probe_table = tmp_path / "probe.tsv"
    probe_table.write_text(
        "id\tFlank\nprobe1\tCGGCAA[Y]GACGCATTCG\n",
        encoding="utf-8",
    )

    offsets, fasta = realign_blast.fasta_and_offsets_from_probe_table(probe_table)

    assert fasta.read_text(encoding="utf-8") == ">probe1\nCGGCAACGACGCATTCG\n"
    assert offsets.to_dict("records") == [
        {"id": "probe1", "offset_fwd": 6, "offset_rev": 10, "alleles": "C/T"}
    ]


def test_fasta_and_offsets_from_probe_table_keeps_slash_marker_and_flank_iupac(
    tmp_path: Path,
) -> None:
    probe_table = tmp_path / "probe.tsv"
    probe_table.write_text("id\tFlank\nprobe1\tARCGT[A/G]TY\n", encoding="utf-8")

    offsets, fasta = realign_blast.fasta_and_offsets_from_probe_table(probe_table)

    assert fasta.read_text(encoding="utf-8") == ">probe1\nAACGTATC\n"
    assert offsets.to_dict("records") == [
        {"id": "probe1", "offset_fwd": 5, "offset_rev": 2, "alleles": "A/G"}
    ]


def test_fasta_and_offsets_from_probe_table_expands_multiple_markers(tmp_path: Path) -> None:
    probe_table = tmp_path / "probe.tsv"
    probe_table.write_text("id\tFlank\nprobe1\tAC[A/G]TG[A/C]AAA\n", encoding="utf-8")

    offsets, fasta = realign_blast.fasta_and_offsets_from_probe_table(probe_table)

    assert fasta.read_text(encoding="utf-8") == ">probe1\nACATGAAAA\n"
    assert offsets.to_dict("records") == [
        {"id": "probe1", "offset_fwd": 2, "offset_rev": 6, "alleles": "A/G"},
        {"id": "probe1", "offset_fwd": 5, "offset_rev": 3, "alleles": "A/C"},
    ]


def test_fasta_and_offsets_from_probe_table_reports_context_for_invalid_marker(
    tmp_path: Path,
) -> None:
    probe_table = tmp_path / "probe.tsv"
    probe_table.write_text("id\tFlank\nprobe1\tACGT[Z]TGCA\n", encoding="utf-8")

    with pytest.raises(ValueError, match="probe1.*ACGT\\[Z\\]TGCA"):
        realign_blast.fasta_and_offsets_from_probe_table(probe_table)


def test_load_id_chrom_map_allows_multiple_target_chroms(tmp_path: Path) -> None:
    id_chrom_map = tmp_path / "id_chrom.tsv"
    id_chrom_map.write_text("probe1\tchr1\nprobe1\tchr2\nprobe2\tchr3\n", encoding="utf-8")

    assert realign_blast.load_id_chrom_map(id_chrom_map) == {
        "probe1": {"chr1", "chr2"},
        "probe2": {"chr3"},
    }


def test_build_id_mapping_strictly_filters_by_id_chrom_map() -> None:
    alignments = pd.DataFrame(
        [
            _alignment_row(chrom="chr1", bitscore=200),
            _alignment_row(chrom="chr2", bitscore=100),
        ]
    )

    result = realign_blast.build_id_mapping(
        alignments,
        _offsets(),
        id_chrom_map={"probe1": {"chr2"}},
        max_gap_opens=0,
    )

    assert result[["id", "chrom", "pos"]].to_dict("records") == [
        {"id": "probe1", "chrom": "chr2", "pos": 105}
    ]


def test_build_id_mapping_drops_id_without_allowed_target_chrom() -> None:
    alignments = pd.DataFrame([_alignment_row(chrom="chr1")])

    result = realign_blast.build_id_mapping(
        alignments,
        _offsets(),
        id_chrom_map={"probe1": {"chr2"}},
    )

    assert result.empty
    assert {"id", "new_id", "pos", "alleles"}.issubset(result.columns)


def test_build_id_mapping_outputs_multiple_rows_for_multiple_marker_offsets() -> None:
    alignments = pd.DataFrame([_alignment_row()])
    offsets = pd.DataFrame(
        {
            "id": ["probe1", "probe1"],
            "offset_fwd": [2, 5],
            "offset_rev": [7, 4],
            "alleles": ["A/G", "A/C"],
        }
    )

    result = realign_blast.build_id_mapping(
        alignments,
        offsets,
        max_gap_opens=0,
    )

    assert result[["id", "chrom", "pos", "alleles"]].to_dict("records") == [
        {"id": "probe1", "chrom": "chr1", "pos": 102, "alleles": "A/G"},
        {"id": "probe1", "chrom": "chr1", "pos": 105, "alleles": "A/C"},
    ]


def test_build_id_mapping_uses_btop_to_correct_gap_before_offset() -> None:
    alignments = pd.DataFrame(
        [
            _alignment_row(
                gap_opens=1,
                btop="3-A7",
                send=110,
                align_len=11,
            )
        ]
    )

    result = realign_blast.build_id_mapping(
        alignments,
        _offsets(offset=5),
        max_gap_opens=2,
    )

    assert result[["id", "chrom", "pos"]].to_dict("records") == [
        {"id": "probe1", "chrom": "chr1", "pos": 106}
    ]


def test_build_id_mapping_drops_hit_when_target_offset_is_on_subject_gap() -> None:
    alignments = pd.DataFrame(
        [
            _alignment_row(
                gap_opens=1,
                btop="5A-4",
                send=108,
                align_len=10,
            )
        ]
    )

    result = realign_blast.build_id_mapping(
        alignments,
        _offsets(offset=5),
        max_gap_opens=2,
    )

    assert result.empty


def test_build_id_mapping_requires_btop_when_gaps_are_allowed() -> None:
    alignments = pd.DataFrame([_alignment_row(gap_opens=1, btop=pd.NA)])

    with pytest.raises(ValueError, match="BTOP.*--force"):
        realign_blast.build_id_mapping(
            alignments,
            _offsets(),
            max_gap_opens=2,
        )


def test_build_id_mapping_with_btop_handles_negative_strand() -> None:
    alignments = pd.DataFrame([_alignment_row(sstart=200, send=191)])

    result = realign_blast.build_id_mapping(
        alignments,
        _offsets(offset=3),
        max_gap_opens=0,
    )

    assert result[["id", "chrom", "pos"]].to_dict("records") == [
        {"id": "probe1", "chrom": "chr1", "pos": 197}
    ]
