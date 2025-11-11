from pathlib import Path
from unittest import mock

import pandas as pd

from liftover import liftover_by_vcf


def test_infer_probe_name_handles_compressed_suffixes(tmp_path):
    gz_vcf = tmp_path / "variants.vcf.gz"
    gz_vcf.touch()
    name = liftover_by_vcf.infer_probe_name(gz_vcf)
    assert name == "variants"

    bcf_vcf = tmp_path / "sample.bcf"
    bcf_vcf.touch()
    assert liftover_by_vcf.infer_probe_name(bcf_vcf) == "sample"


def test_build_output_paths_uses_probe_name(tmp_path):
    vcf = tmp_path / "input.vcf"
    vcf.touch()
    outdir = tmp_path / "out"
    outdir.mkdir()

    outputs = liftover_by_vcf.build_output_paths(vcf, outdir)
    assert outputs.probe_name == "input"
    assert outputs.probe_bed == outdir / "input.bed"
    assert outputs.raw_bed == outdir / "raw.input.bed"


def test_liftover_vcf_pipeline_invokes_external_tools(tmp_path):
    vcf = tmp_path / "input.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\tid\tA\tT\t.\t.\t.\n",
        encoding="utf-8",
    )
    chain = tmp_path / "hg19ToHg38.over.chain.gz"
    chain.touch()
    ref_fa = tmp_path / "ref.fa"
    ref_fa.write_text(">chr1\nAAAAAAAAAA\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t10\t4\t10\t11\n", encoding="utf-8")
    query_fa = tmp_path / "query.fa"
    query_fa.write_text(">chr1\nAAAAAAAAAA\n", encoding="utf-8")
    (tmp_path / "query.fa.fai").write_text("chr1\t10\t4\t10\t11\n", encoding="utf-8")
    outdir = tmp_path / "lift"

    mock_df = pd.DataFrame({"chrom": ["chr1"], "pos": [10]})
    mock_df["start"] = mock_df["pos"] - 1
    mock_df["id"] = mock_df["chrom"] + "_" + mock_df["pos"].astype(str)

    with mock.patch("liftover.liftover_by_vcf.delegator.run") as mock_run, mock.patch(
        "liftover.liftover_by_vcf.read_lifted_positions", return_value=mock_df
    ):
        liftover_by_vcf.liftover_vcf(vcf, chain, ref_fa, query_fa, outdir, force=True)

    assert mock_run.call_count == 3  # liftvcf, bedtools sort, bedtools slop+merge
    assert (outdir / "input.bed").exists()
    assert (outdir / "input.id").exists()
    assert (outdir / "input.snpcalling.bed").exists()
