import importlib.util
from pathlib import Path

import pandas as pd
import pytest


MODULE_PATH = Path(__file__).resolve().parents[1] / "split-genome.py"


@pytest.fixture(scope="module")
def split_genome_module():
    spec = importlib.util.spec_from_file_location("split_genome", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_parse_gff_auto_selects_priority_feature(split_genome_module, tmp_path):
    gff_content = "\n".join(
        [
            "chr1\t.\tCDS\t5\t10\t.\t+\t.\tID=cds1",
            "chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1",
            "chr1\t.\tmRNA\t200\t300\t.\t+\t.\tID=mrna1",
            "chr2\t.\tgene\t50\t150\t.\t-\t.\tID=gene2",
        ]
    )
    gff_path = tmp_path / "test.gff"
    gff_path.write_text(gff_content)

    result = split_genome_module.parse_gff_with_pandas(gff_path)

    assert not result.empty
    assert result["feature"].eq("gene").all()
    assert set(result["attributes"]) == {"ID=gene1", "ID=gene2"}


def test_select_split_candidates_returns_best_gap(split_genome_module):
    gap_df = pd.DataFrame(
        [
            {"Chromosome": "chr1", "gap_start": 100, "gap_end": 200, "gap_size": 101},
            {"Chromosome": "chr1", "gap_start": 210, "gap_end": 260, "gap_size": 51},
            {"Chromosome": "chr2", "gap_start": 50, "gap_end": 120, "gap_size": 71},
        ]
    )
    chrom_bounds = pd.DataFrame(
        [
            {
                "Chromosome": "chr1",
                "chrom_size": 800,
                "split_start": 90,
                "split_end": 600,
            },
            {
                "Chromosome": "chr2",
                "chrom_size": 400,
                "split_start": 40,
                "split_end": 300,
            },
        ]
    )

    best = split_genome_module.select_split_candidates(chrom_bounds, gap_df, 50)

    assert len(best) == 2
    chr1_row = best.loc[best["Chromosome"] == "chr1"].iloc[0]
    assert chr1_row["gap_start"] == 100
    assert chr1_row["gap_end"] == 200


def test_select_split_candidates_raises_for_small_gap(split_genome_module):
    gap_df = pd.DataFrame(
        [
            {"Chromosome": "chr1", "gap_start": 100, "gap_end": 120, "gap_size": 21},
        ]
    )
    chrom_bounds = pd.DataFrame(
        [
            {
                "Chromosome": "chr1",
                "chrom_size": 500,
                "split_start": 50,
                "split_end": 400,
            }
        ]
    )

    with pytest.raises(ValueError) as excinfo:
        split_genome_module.select_split_candidates(chrom_bounds, gap_df, 100)

    assert "小于最小阈值" in str(excinfo.value)


def test_generate_split_genome_outputs_expected_sequences(
    split_genome_module, tmp_path
):
    fasta_content = ">chr1\nAAAACCCC\n>chr2\nGGGGG\n"
    fasta_path = tmp_path / "genome.fa"
    fasta_path.write_text(fasta_content)

    bed_df = pd.DataFrame(
        [
            {"Chromosome": "chr1", "start": 0, "end": 4, "id": "chr1a"},
            {"Chromosome": "chr1", "start": 4, "end": 8, "id": "chr1b"},
            {"Chromosome": "chr2", "start": 0, "end": 5, "id": "chr2"},
        ]
    )
    out_path = tmp_path / "genome.split.fa"

    split_genome_module.generate_split_genome(
        fasta_path, bed_df, out_path, show_progress=False
    )

    assert out_path.read_text().strip().splitlines() == [
        ">chr1a",
        "AAAA",
        ">chr1b",
        "CCCC",
        ">chr2",
        "GGGGG",
    ]
