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


def test_build_id_mapping_passes_through_id_missing_from_id_chrom_map() -> None:
    alignments = pd.DataFrame(
        [
            _alignment_row(probe_id="probe1", chrom="chr2"),
            _alignment_row(probe_id="probe2", chrom="chr9"),
        ]
    )
    offsets = pd.DataFrame(
        {
            "id": ["probe1", "probe2"],
            "offset_fwd": [5, 5],
            "offset_rev": [4, 4],
            "alleles": ["A/G", "A/T"],
        }
    )

    result = realign_blast.build_id_mapping(
        alignments,
        offsets,
        id_chrom_map={"probe1": {"chr2"}},
        max_gap_opens=0,
    )

    assert sorted(result["id"].tolist()) == ["probe1", "probe2"]
    assert result[result["id"] == "probe2"]["chrom"].tolist() == ["chr9"]


def test_build_id_mapping_still_filters_mapped_id_when_other_ids_missing() -> None:
    alignments = pd.DataFrame(
        [
            _alignment_row(probe_id="probe1", chrom="chr1"),
            _alignment_row(probe_id="probe2", chrom="chr9"),
        ]
    )
    offsets = pd.DataFrame(
        {
            "id": ["probe1", "probe2"],
            "offset_fwd": [5, 5],
            "offset_rev": [4, 4],
            "alleles": ["A/G", "A/T"],
        }
    )

    result = realign_blast.build_id_mapping(
        alignments,
        offsets,
        id_chrom_map={"probe1": {"chr2"}},
        max_gap_opens=0,
    )

    # probe1 被严格过滤掉（chr1 不在允许集合），probe2 未登记则放行
    assert result["id"].tolist() == ["probe2"]


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


def test_fasta_and_offsets_from_probe_table_preserves_n_in_literal_flank(
    tmp_path: Path,
) -> None:
    probe_table = tmp_path / "probe.tsv"
    probe_table.write_text(
        "id\tFlank\nprobe1\tNNAC[A/G]TGNN\n",
        encoding="utf-8",
    )

    offsets, fasta = realign_blast.fasta_and_offsets_from_probe_table(probe_table)

    assert fasta.read_text(encoding="utf-8") == ">probe1\nNNACATGNN\n"
    assert offsets.to_dict("records") == [
        {"id": "probe1", "offset_fwd": 4, "offset_rev": 4, "alleles": "A/G"}
    ]


def test_count_ns_in_fasta_handles_multiline_sequences(tmp_path: Path) -> None:
    fa = tmp_path / "t.fa"
    fa.write_text(">a\nACNN\nNNAT\n>b\nACGT\n>c\nNnNn\n", encoding="utf-8")

    assert realign_blast.count_ns_in_fasta(fa) == {"a": 4, "b": 0, "c": 4}


def test_build_id_mapping_uses_informative_len_for_match_ratio() -> None:
    # query_len=10, align_len=8, n_count=3 → informative_len=7 → ratio=8/7 > 0.9 → 通过
    alignments = pd.DataFrame([_alignment_row(align_len=8, send=107)])

    result = realign_blast.build_id_mapping(
        alignments,
        _offsets(),
        max_gap_opens=0,
        match_ratio_cutoff=0.9,
        n_counts={"probe1": 3},
    )

    assert result["id"].tolist() == ["probe1"]


def test_build_id_mapping_still_rejects_low_ratio_when_no_n() -> None:
    # 无 N：align_len=5, query_len=10 → ratio=0.5 → 不过
    alignments = pd.DataFrame([_alignment_row(align_len=5, send=104)])

    result = realign_blast.build_id_mapping(
        alignments,
        _offsets(),
        max_gap_opens=0,
        match_ratio_cutoff=0.9,
    )

    assert result.empty


def test_build_id_mapping_drops_query_with_all_n_sequence() -> None:
    alignments = pd.DataFrame([_alignment_row()])

    result = realign_blast.build_id_mapping(
        alignments,
        _offsets(),
        max_gap_opens=0,
        n_counts={"probe1": 10},  # informative_len ≤ 0
    )

    assert result.empty


def test_build_id_mapping_annotates_rank_and_status_columns() -> None:
    alignments = pd.DataFrame(
        [
            _alignment_row(chrom="chr1", bitscore=200),
            _alignment_row(chrom="chr2", bitscore=100),
        ]
    )

    result = realign_blast.build_id_mapping(
        alignments,
        _offsets(),
        max_gap_opens=0,
        max_hits=2,
    )

    assert result["rank"].tolist() == [1, 2]
    assert result["chrom"].tolist() == ["chr1", "chr2"]
    assert result["chr_map_status"].tolist() == ["n/a", "n/a"]
    assert result["id_chrom_status"].tolist() == ["n/a", "n/a"]
    assert result["selection_reason"].str.contains("bitscore").all()
    assert result["selection_reason"].iloc[0].startswith("rank #1")


def test_build_id_mapping_status_is_strict_when_id_is_registered() -> None:
    alignments = pd.DataFrame([_alignment_row(chrom="chr1")])

    result = realign_blast.build_id_mapping(
        alignments,
        _offsets(),
        max_gap_opens=0,
        id_chrom_map={"probe1": {"chr1"}},
    )

    assert result["id_chrom_status"].tolist() == ["strict"]
    assert "id_chrom_map" in result["selection_reason"].iloc[0]


def test_build_id_mapping_status_is_fallback_when_id_missing_from_chrom_map() -> None:
    alignments = pd.DataFrame([_alignment_row(chrom="chr1")])

    result = realign_blast.build_id_mapping(
        alignments,
        _offsets(),
        max_gap_opens=0,
        id_chrom_map={"probe2": {"chr1"}},
    )

    assert result["id_chrom_status"].tolist() == ["fallback"]


def test_build_id_mapping_chr_map_status_reflects_match_vs_fallback() -> None:
    alignments = pd.DataFrame(
        [
            _alignment_row(chrom="chrTarget", bitscore=200),
            _alignment_row(chrom="chrOther", bitscore=150),
        ]
    )
    source_chroms = pd.DataFrame({"id": ["probe1"], "source_chrom": ["chrSrc"]})

    result = realign_blast.build_id_mapping(
        alignments,
        _offsets(),
        max_gap_opens=0,
        max_hits=2,
        chr_map={"chrSrc": {"chrTarget"}},
        source_chroms=source_chroms,
    )

    # chr_matched=True 排在前；对应 matched；不匹配的 chrOther 记为 fallback
    assert result["chrom"].tolist() == ["chrTarget", "chrOther"]
    assert result["chr_map_status"].tolist() == ["matched", "fallback"]


def test_write_selection_report_contains_header_and_reasons(tmp_path: Path) -> None:
    alignments = pd.DataFrame([_alignment_row()])
    mapping_df = realign_blast.build_id_mapping(
        alignments,
        _offsets(),
        max_gap_opens=0,
    )

    report_path = realign_blast.write_selection_report(mapping_df, tmp_path / "out")

    assert report_path.exists()
    content = report_path.read_text(encoding="utf-8").strip().splitlines()
    header = content[0].split("\t")
    assert header[0] == "id"
    assert "selection_reason" in header
    assert "match_ratio" in header
    assert len(content) == 2  # 表头 + 1 条记录


def test_write_selection_report_drops_all_n_a_status_columns(tmp_path: Path) -> None:
    # 未传 chr_map / id_chrom_map：两列状态均为 n/a，应被剔除
    alignments = pd.DataFrame([_alignment_row()])
    mapping_df = realign_blast.build_id_mapping(alignments, _offsets(), max_gap_opens=0)

    report_path = realign_blast.write_selection_report(mapping_df, tmp_path / "out")

    header = report_path.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert "chr_map_status" not in header
    assert "id_chrom_status" not in header


def test_write_selection_report_keeps_status_when_values_vary(tmp_path: Path) -> None:
    # id_chrom_map 下状态为 strict：id_chrom_status 应保留；chr_map_status 仍全 n/a 被剔除
    alignments = pd.DataFrame([_alignment_row(chrom="chr1")])
    mapping_df = realign_blast.build_id_mapping(
        alignments,
        _offsets(),
        max_gap_opens=0,
        id_chrom_map={"probe1": {"chr1"}},
    )

    report_path = realign_blast.write_selection_report(mapping_df, tmp_path / "out")

    header = report_path.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert "id_chrom_status" in header
    assert "chr_map_status" not in header


def test_write_outputs_uses_separate_selection_df_for_report(tmp_path: Path) -> None:
    alignments = pd.DataFrame(
        [
            _alignment_row(chrom="chr1", bitscore=200),
            _alignment_row(chrom="chr2", bitscore=100),
        ]
    )
    full_df = realign_blast.build_id_mapping(
        alignments,
        _offsets(),
        max_gap_opens=0,
        max_hits=2,
    )
    primary_df = full_df[full_df["rank"] <= 1].copy()

    out_prefix = tmp_path / "out"
    realign_blast.write_outputs(primary_df, out_prefix, selection_df=full_df)

    idmap_lines = (tmp_path / "out.idmap.tsv").read_text(encoding="utf-8").strip().splitlines()
    selection_lines = (
        (tmp_path / "out.idmap.selection.tsv").read_text(encoding="utf-8").strip().splitlines()
    )

    # idmap.tsv 不带表头：仅保留 rank=1 的一条
    assert len(idmap_lines) == 1
    # selection.tsv 带表头：表头 + rank=1 与 rank=2 两条
    assert len(selection_lines) == 3
    assert selection_lines[0].split("\t")[0] == "id"


def test_write_outputs_defaults_selection_df_to_primary(tmp_path: Path) -> None:
    alignments = pd.DataFrame([_alignment_row()])
    mapping_df = realign_blast.build_id_mapping(
        alignments,
        _offsets(),
        max_gap_opens=0,
    )

    realign_blast.write_outputs(mapping_df, tmp_path / "out")

    selection_lines = (
        (tmp_path / "out.idmap.selection.tsv").read_text(encoding="utf-8").strip().splitlines()
    )
    assert len(selection_lines) == 2  # 表头 + 1 行
