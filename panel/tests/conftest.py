from pathlib import Path

import pandas as pd
import pytest


@pytest.fixture
def sample_blast_data():
    """创建示例 BLAST 数据"""
    return pd.DataFrame(
        {
            "id": ["seq1", "seq1", "seq2", "seq3"],
            "identity": [98, 95, 100, 92],
            "match_len": [100, 80, 120, 90],
            "mismatch": [2, 4, 0, 7],
            "gapopen": [0, 1, 0, 1],
        }
    )


@pytest.fixture
def temp_blast_dir(tmp_path):
    """创建临时 BLAST 文件目录"""
    blast_dir = tmp_path / "blast_results"
    blast_dir.mkdir()

    sample_data = (
        "seq1\t1e-50\t98\t100\t2\t0\n"
        "seq1\t1e-40\t95\t80\t4\t1\n"
        "seq2\t1e-60\t100\t120\t0\t0\n"
    )

    test_file = blast_dir / "test1.blast"
    test_file.write_text(sample_data)

    return blast_dir


@pytest.fixture
def blast_config():
    """创建测试配置"""
    from panel.blast_res_summary import BlastConfig

    return BlastConfig(
        min_match_length=24,
        max_match_count=1,
        probe_length=120,
        output_all=False,
        show_max_length=False,
        min_identity=90,
        max_mismatch_count=5,
        max_gap_count=2,
    )
