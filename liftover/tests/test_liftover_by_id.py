import gzip
from pathlib import Path
from unittest import mock

import pandas as pd
import pytest
from pyfaidx import Faidx, Fasta

from liftover import liftover_by_id


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_fasta(tmp_path: Path, name: str = "ref.fa", seqs: dict[str, str] | None = None) -> Path:
    """写入一个 pyfaidx 可索引的 FASTA 文件并创建 .fai 索引。"""
    if seqs is None:
        seqs = {"chr1": "ACGTACGTAC"}
    fa_path = tmp_path / name
    with open(fa_path, "w", encoding="utf-8") as fh:
        for chrom, seq in seqs.items():
            fh.write(f">{chrom}\n{seq}\n")
    Faidx(str(fa_path))  # 生成 .fai
    return fa_path


def _write_gzipped_vcf(path: Path, records: list[str]) -> Path:
    """写入一个 gzip 压缩的最小 VCF 文件。"""
    lines = [
        "##fileformat=VCFv4.2\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    ]
    lines.extend(f"{r}\n" for r in records)
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.writelines(lines)
    return path


# ---------------------------------------------------------------------------
# parse_id_to_chrom_pos
# ---------------------------------------------------------------------------

class TestParseIdToChromPos:
    def test_simple_id(self):
        assert liftover_by_id.parse_id_to_chrom_pos("chr1_12345") == ("chr1", 12345)

    def test_multi_underscore_id(self):
        assert liftover_by_id.parse_id_to_chrom_pos("chr1_random_500") == ("chr1_random", 500)

    def test_no_underscore_raises(self):
        with pytest.raises(ValueError):
            liftover_by_id.parse_id_to_chrom_pos("nounderscore")

    def test_non_numeric_pos_raises(self):
        with pytest.raises(ValueError):
            liftover_by_id.parse_id_to_chrom_pos("chr1_abc")


# ---------------------------------------------------------------------------
# fetch_ref_nucleotide
# ---------------------------------------------------------------------------

class TestFetchRefNucleotide:
    def test_returns_correct_base(self, tmp_path: Path):
        fa_path = _write_fasta(tmp_path, seqs={"chr1": "ACGTACGTAC"})
        ref_fa = Fasta(str(fa_path))
        assert liftover_by_id.fetch_ref_nucleotide(ref_fa, "chr1", 1) == "A"
        assert liftover_by_id.fetch_ref_nucleotide(ref_fa, "chr1", 4) == "T"

    def test_missing_chrom_returns_n(self, tmp_path: Path):
        fa_path = _write_fasta(tmp_path, seqs={"chr1": "ACGT"})
        ref_fa = Fasta(str(fa_path))
        assert liftover_by_id.fetch_ref_nucleotide(ref_fa, "chrZ", 1) == "N"


# ---------------------------------------------------------------------------
# read_liftover_result
# ---------------------------------------------------------------------------

class TestReadLiftoverResult:
    def test_reads_vcf_gz_and_deduplicates(self, tmp_path: Path):
        vcf_gz = tmp_path / "lifted.vcf.gz"
        _write_gzipped_vcf(
            vcf_gz,
            [
                "chr1\t100\tinput_chr1_1\tA\tT\t.\t.\t.",
                "chr1\t200\tinput_chr1_2\tC\tG\t.\t.\t.",
                "chr1\t100\tinput_chr1_1\tA\tT\t.\t.\t.",  # duplicate
            ],
        )
        df = liftover_by_id.read_liftover_result(vcf_gz)
        assert len(df) == 2
        assert "id" in df.columns
        assert "start" in df.columns
        assert "pos_id" in df.columns
        assert df.iloc[0]["id"] == "input_chr1_1"
        assert df.iloc[0]["start"] == 99
        assert df.iloc[0]["pos_id"] == "chr1_100"


# ---------------------------------------------------------------------------
# load_chrom_sizes_from_fai
# ---------------------------------------------------------------------------

class TestLoadChromSizesFromFai:
    def test_reads_fai(self, tmp_path: Path):
        fai = tmp_path / "genome.fa.fai"
        fai.write_text(
            "chr1\t1000\t5\t80\t81\n"
            "chr2\t500\t1100\t80\t81\n",
            encoding="utf-8",
        )
        df = liftover_by_id.load_chrom_sizes_from_fai(fai)
        assert list(df.columns) == ["chrom", "chrom_size"]
        assert len(df) == 2
        assert df.iloc[0]["chrom"] == "chr1"
        assert df.iloc[0]["chrom_size"] == 1000


# ---------------------------------------------------------------------------
# sort_bed_by_fai
# ---------------------------------------------------------------------------

class TestSortBedByFai:
    def test_sorts_by_chrom_order_and_start(self):
        bed_df = pd.DataFrame({
            "chrom": ["chr2", "chr1", "chr1"],
            "start": [50, 200, 100],
            "end": [51, 201, 101],
        })
        chrom_df = pd.DataFrame({"chrom": ["chr1", "chr2"], "chrom_size": [1000, 500]})
        result = liftover_by_id.sort_bed_by_fai(bed_df, chrom_df)
        assert list(result["chrom"]) == ["chr1", "chr1", "chr2"]
        assert list(result["start"]) == [100, 200, 50]

    def test_empty_dataframe_returns_empty(self):
        bed_df = pd.DataFrame(columns=["chrom", "start", "end"])
        chrom_df = pd.DataFrame({"chrom": ["chr1"], "chrom_size": [1000]})
        result = liftover_by_id.sort_bed_by_fai(bed_df, chrom_df)
        assert len(result) == 0


# ---------------------------------------------------------------------------
# slop_and_merge
# ---------------------------------------------------------------------------

class TestSlopAndMerge:
    def test_two_distant_snps_stay_separate(self):
        bed_df = pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "start": [99, 999],
            "end": [100, 1000],
        })
        chrom_sizes = {"chr1": 5000}
        result = liftover_by_id.slop_and_merge(bed_df, chrom_sizes, flank=100)
        assert len(result) == 2

    def test_two_close_snps_merge(self):
        bed_df = pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "start": [99, 149],
            "end": [100, 150],
        })
        chrom_sizes = {"chr1": 5000}
        result = liftover_by_id.slop_and_merge(bed_df, chrom_sizes, flank=100)
        assert len(result) == 1
        assert result.iloc[0]["start"] == 0  # 99 - 100 = -1 → clamped to 0
        assert result.iloc[0]["end"] == 250  # 150 + 100

    def test_clamp_to_zero_near_start(self):
        bed_df = pd.DataFrame({"chrom": ["chr1"], "start": [9], "end": [10]})
        chrom_sizes = {"chr1": 5000}
        result = liftover_by_id.slop_and_merge(bed_df, chrom_sizes, flank=100)
        assert result.iloc[0]["start"] == 0

    def test_clamp_to_chrom_size_near_end(self):
        bed_df = pd.DataFrame({"chrom": ["chr1"], "start": [4949], "end": [4950]})
        chrom_sizes = {"chr1": 5000}
        result = liftover_by_id.slop_and_merge(bed_df, chrom_sizes, flank=100)
        assert result.iloc[0]["end"] == 5000

    def test_empty_input_returns_empty(self):
        bed_df = pd.DataFrame(columns=["chrom", "start", "end"])
        result = liftover_by_id.slop_and_merge(bed_df, {}, flank=100)
        assert len(result) == 0


# ---------------------------------------------------------------------------
# make_id_vcf
# ---------------------------------------------------------------------------

class TestMakeIdVcf:
    def test_creates_vcf_from_id_file(self, tmp_path: Path):
        fa_path = _write_fasta(tmp_path, seqs={"chr1": "ACGTACGTAC"})
        id_file = tmp_path / "snps.id"
        id_file.write_text("chr1_1\nchr1_4\n", encoding="utf-8")

        vcf_path = liftover_by_id.make_id_vcf(id_file, fa_path, force=False)
        assert vcf_path.exists()
        content = vcf_path.read_text(encoding="utf-8")
        assert "##fileformat=VCFv4.2" in content
        assert "chr1\t1\tchr1_1\tA\t" in content
        assert "chr1\t4\tchr1_4\tT\t" in content

    def test_skip_when_exists_and_no_force(self, tmp_path: Path):
        fa_path = _write_fasta(tmp_path, seqs={"chr1": "ACGT"})
        id_file = tmp_path / "snps.id"
        id_file.write_text("chr1_1\n", encoding="utf-8")

        # 先创建一次
        vcf_path = liftover_by_id.make_id_vcf(id_file, fa_path, force=False)
        mtime_first = vcf_path.stat().st_mtime

        # 再次调用，force=False 应跳过
        vcf_path2 = liftover_by_id.make_id_vcf(id_file, fa_path, force=False)
        assert vcf_path2.stat().st_mtime == mtime_first


# ---------------------------------------------------------------------------
# Integration: main CLI (mock subprocess.run for transanno)
# ---------------------------------------------------------------------------

def test_main_pipeline_invokes_subprocess_and_creates_outputs(tmp_path: Path):
    """验证 main 调用 subprocess.run 并生成预期输出文件。"""
    id_file = tmp_path / "test.id"
    id_file.write_text("chr1_10\n", encoding="utf-8")

    chain = tmp_path / "hg19ToHg38.over.chain.gz"
    chain.touch()

    ref_fa = _write_fasta(tmp_path, name="ref.fa", seqs={"chr1": "AAAAAAAAAAAA"})
    query_fa = _write_fasta(tmp_path, name="query.fa", seqs={"chr1": "AAAAAAAAAAAA"})

    outdir = tmp_path / "lift"

    mock_df = pd.DataFrame({"chrom": ["chr1"], "pos": [10], "id": ["chr1_10"]})
    mock_df["start"] = mock_df["pos"] - 1
    mock_df["pos_id"] = mock_df["chrom"] + "_" + mock_df["pos"].astype(str)

    mock_proc = mock.MagicMock(returncode=0, stderr="")
    with (
        mock.patch("liftover.liftover_by_id.subprocess.run", return_value=mock_proc) as mock_run,
        mock.patch(
            "liftover.liftover_by_id.read_liftover_result", return_value=mock_df,
        ),
    ):
        liftover_by_id.main(
            id_file=id_file,
            chain=chain,
            ref_fa=ref_fa,
            query_fa=query_fa,
            outdir=outdir,
            force=True,
        )

    assert mock_run.call_count >= 1
    assert (outdir / "test.id").exists()
    assert (outdir / "test.bed").exists()
    assert (outdir / "test.pos.tsv").exists()
    assert (outdir / "test.snpcalling.bed").exists()
    assert (outdir / "test.pos.tsv").read_text(encoding="utf-8").strip() == "chr1\t10\tchr1_10"
