from io import StringIO

import pandas as pd

from gtf.id_map_from_gffcompare import (
    generate_class_code_df,
    load_gffcompare_tmap,
    process_overlap_mapping,
)


def test_generate_class_code_df():
    df = generate_class_code_df()
    assert list(df["class_code"]) == ["=", "k", "c"]
    assert list(df["rank"]) == [0, 1, 2]


def test_process_overlap_mapping():
    gff_data = """qry_gene_id\tclass_code\tref_gene_id
q1\t=\tr1
q2\tk\tr2
q3\tc\tr3
q4\to\tr4
"""
    gff_df = pd.read_csv(StringIO(gff_data), sep="\t")

    group_data = """group_id\tiwgsc_v2.1
g1\tr1
g2\tr2
g3\trX
"""
    group_df = pd.read_csv(StringIO(group_data), sep="\t")

    result_df = process_overlap_mapping(gff_df, group_df, name="mysample")

    expected = pd.DataFrame(
        {
            "group_id": ["g1", "g2"],
            "mysample": ["q1", "q2"],
        }
    )

    pd.testing.assert_frame_equal(result_df.reset_index(drop=True), expected)


def test_load_gffcompare_tmap(tmp_path):
    content = "qry_gene_id\tX\tclass_code\tref_gene_id\nq1\t.\t=\tr1\n"
    file = tmp_path / "test.tmap"
    file.write_text(content)

    df = load_gffcompare_tmap(file)
    assert list(df.columns) == ["qry_gene_id", "class_code", "ref_gene_id"]
    assert df.iloc[0]["ref_gene_id"] == "r1"
    file = tmp_path / "test.tmap"
    file.write_text(content)

    df = load_gffcompare_tmap(file)
    assert list(df.columns) == ["qry_gene_id", "class_code", "ref_gene_id"]
    assert df.iloc[0]["ref_gene_id"] == "r1"
