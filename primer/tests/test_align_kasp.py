#!/usr/bin/env python3
import importlib.util
from pathlib import Path

import pandas as pd
import pytest


@pytest.fixture(scope="session")
def align_kasp_module():
    script_path = Path(__file__).parent.parent / "align-kasp.py"
    spec = importlib.util.spec_from_file_location("primer_align_kasp", script_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load module from {script_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_iupac_to_atgc_basic_and_error(align_kasp_module):
    assert align_kasp_module.IUPAC_to_ATGC("ACGT") == "ACGT"
    assert align_kasp_module.IUPAC_to_ATGC("RY") == "AC"  # R->A, Y->C
    with pytest.raises(ValueError, match="不支持的IUPAC码"):
        align_kasp_module.IUPAC_to_ATGC("Z")


def test_ssr_table_to_fa_happy_path(tmp_path, align_kasp_module):
    df = pd.DataFrame(
        {
            "name": ["SSR1", "SSR2"],
            "left": ["ATCG", "GCTA"],
            "right": ["CGAT", "TAGC"],
        }
    )
    left_fa, right_fa = align_kasp_module.ssr_table_to_fa(df, tmp_path, prefix="t")
    assert left_fa.read_text() == ">SSR1\nATCG\n>SSR2\nGCTA\n"
    assert right_fa.read_text() == ">SSR1\nCGAT\n>SSR2\nTAGC\n"


def test_ssr_table_to_fa_validates_inputs(tmp_path, align_kasp_module):
    df_missing = pd.DataFrame({"name": ["SSR1"], "left": ["ATCG"]})
    with pytest.raises(ValueError, match="缺少必需的列"):
        align_kasp_module.ssr_table_to_fa(df_missing, tmp_path)

    with pytest.raises(FileNotFoundError, match="输出目录不存在"):
        align_kasp_module.ssr_table_to_fa(
            pd.DataFrame({"name": ["SSR1"], "left": ["ATCG"], "right": ["CGAT"]}),
            tmp_path / "nope",
        )


def test_kasp_table_to_fa_happy_path(tmp_path, align_kasp_module):
    df = pd.DataFrame(
        {
            "name": ["K1"],
            "fam": ["ATCG"],
            "vic": ["GCTA"],
            "common": ["CGAT"],
        }
    )
    fam_fa, vic_fa, common_fa = align_kasp_module.kasp_table_to_fa(
        df, tmp_path, prefix="t"
    )
    assert fam_fa.read_text() == ">K1\nATCG\n"
    assert vic_fa.read_text() == ">K1\nGCTA\n"
    assert common_fa.read_text() == ">K1\nCGAT\n"


def test_best_match_pos_selects_best(align_kasp_module):
    df = pd.DataFrame(
        [
            {
                "name": "a",
                "left_mismatch": 1,
                "right_mismatch": 0,
                "left_gap": 0,
                "right_gap": 0,
                "left_match_length": 10,
                "right_match_length": 10,
                "extra": "keep1",
            },
            {
                "name": "a",
                "left_mismatch": 0,
                "right_mismatch": 0,
                "left_gap": 1,
                "right_gap": 0,
                "left_match_length": 11,
                "right_match_length": 10,
                "extra": "keep2",
            },
            {
                "name": "a",
                "left_mismatch": 0,
                "right_mismatch": 0,
                "left_gap": 0,
                "right_gap": 0,
                "left_match_length": 11,
                "right_match_length": 10,
                "extra": "best",
            },
        ]
    )
    best = align_kasp_module.best_match_pos(df)
    assert len(best) == 1
    assert best.iloc[0]["extra"] == "best"


def test_align_ssr_seq_skips_when_output_exists(tmp_path, align_kasp_module, monkeypatch):
    fa = tmp_path / "x.fa"
    fa.write_text(">x\nATCG\n")
    blast_db = tmp_path / "db"
    blast_db.write_text("fake")

    blast_out = fa.with_suffix(".blasttab.tsv")
    blast_out.write_text("already\n")

    def _fail_run(*args, **kwargs):
        raise AssertionError("blast runner should not be called when blast output exists")

    if align_kasp_module.delegator is not None:
        monkeypatch.setattr(align_kasp_module.delegator, "run", _fail_run)
    else:
        monkeypatch.setattr(align_kasp_module.subprocess, "run", _fail_run)

    out = align_kasp_module.align_ssr_seq(fa, blast_db, force=False)
    assert out == blast_out


def test_blast_processor_load_blast_df_and_process(tmp_path, align_kasp_module):
    processor = align_kasp_module.BlastProcessor(
        config=align_kasp_module.BlastConfig(max_distance=1000, blast_columns=[0, 1, 8, 9])
    )

    left = tmp_path / "left.out"
    right = tmp_path / "right.out"

    # qseqid sseqid ... sstart send ...
    left.write_text("p1\tchr1\t99\t10\t0\t0\t1\t10\t100\t109\t1e-10\t50\n")
    right.write_text("p1\tchr1\t99\t10\t0\t0\t1\t10\t160\t151\t1e-10\t50\n")  # reverse => '-'

    df = processor.process_blast_results(left, right)
    assert not df.empty
    row = df.iloc[0]
    assert row["chrom"] == "chr1"
    assert row["left_strand"] != row["right_strand"]
    assert row["ssr_length"] > 0


def test_main_generates_pos_table(tmp_path, align_kasp_module, monkeypatch):
    kasp_in = tmp_path / "kasp.tsv"
    kasp_in.write_text("m1\tAAAAAAAAAA\tCCCCCCCCCC\tGGGGGGGGGG\n")
    out_dir = tmp_path / "out"
    blast_db = tmp_path / "blast_db"
    blast_db.write_text("fake-db")

    def fake_align_ssr_seq(ssr_fa: Path, _blast_db: Path, threads: int = 1, force: bool = False):
        blast_out = ssr_fa.with_suffix(".blasttab.tsv")
        if ".fam." in ssr_fa.name:
            blast_out.write_text(
                "m1\tchr1\t100\t10\t0\t0\t1\t10\t100\t109\t1e-10\t50\n"
            )
        elif ".vic." in ssr_fa.name:
            blast_out.write_text(
                "m1\tchr1\t100\t10\t0\t0\t1\t10\t105\t109\t1e-10\t50\n"
            )
        else:
            # common on reverse strand to satisfy strand filter
            blast_out.write_text(
                "m1\tchr1\t100\t10\t0\t0\t1\t10\t160\t151\t1e-10\t50\n"
            )
        return blast_out

    monkeypatch.setattr(align_kasp_module, "align_ssr_seq", fake_align_ssr_seq)

    align_kasp_module.main(kasp_in, blast_db, out_dir, threads=1, force=True)

    pos_out = kasp_in.with_suffix(".pos.tsv")
    assert pos_out.exists()
    result = pd.read_table(pos_out)
    assert {"name", "chrom", "snp_pos", "primer_span"}.issubset(result.columns)
    assert result.loc[0, "chrom"] == "chr1"
    assert int(result.loc[0, "snp_pos"]) == 109


def test_main_respects_max_mismatch_filter(tmp_path, align_kasp_module, monkeypatch):
    kasp_in = tmp_path / "kasp.tsv"
    kasp_in.write_text("m1\tAAAAAAAAAA\tCCCCCCCCCC\tGGGGGGGGGG\n")
    out_dir = tmp_path / "out"
    blast_db = tmp_path / "blast_db"
    blast_db.write_text("fake-db")

    def fake_align_ssr_seq(ssr_fa: Path, _blast_db: Path, threads: int = 1, force: bool = False):
        blast_out = ssr_fa.with_suffix(".blasttab.tsv")
        if ".fam." in ssr_fa.name:
            # fam has too many mismatches
            blast_out.write_text(
                "m1\tchr1\t100\t10\t9\t0\t1\t10\t100\t109\t1e-10\t50\n"
            )
        elif ".vic." in ssr_fa.name:
            blast_out.write_text(
                "m1\tchr1\t100\t10\t0\t0\t1\t10\t105\t109\t1e-10\t50\n"
            )
        else:
            blast_out.write_text(
                "m1\tchr1\t100\t10\t0\t0\t1\t10\t160\t151\t1e-10\t50\n"
            )
        return blast_out

    monkeypatch.setattr(align_kasp_module, "align_ssr_seq", fake_align_ssr_seq)

    align_kasp_module.main(kasp_in, blast_db, out_dir, threads=1, force=True, max_mismatch=3)

    pos_out = kasp_in.with_suffix(".pos.tsv")
    result = pd.read_table(pos_out)
    assert pd.isna(result.loc[0, "chrom"])
