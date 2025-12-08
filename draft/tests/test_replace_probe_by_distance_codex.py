import sys
import importlib.util
from pathlib import Path
import pandas as pd
import pytest
from typer.testing import CliRunner

# Load the module dynamically because of the hyphens in the filename
module_path = Path(__file__).parents[1] / "replace-probe-by-distance-codex.py"
spec = importlib.util.spec_from_file_location("replace_probe_module", module_path)
replace_module = importlib.util.module_from_spec(spec)
sys.modules["replace_probe_module"] = replace_module
spec.loader.exec_module(replace_module)

from replace_probe_module import main

import typer
app = typer.Typer()
app.command()(main)

runner = CliRunner()

@pytest.fixture
def replace_data():
    return pd.DataFrame({
        "chrom": ["1", "1"],
        "pos": [1000, 5000],
        "id": ["rs1", "rs2"],
        "target_id": ["t1", "t2"]
    })

@pytest.fixture
def candidate_data():
    return pd.DataFrame({
        "chrom": ["1", "1", "1", "1", "2"],
        "pos": [1010, 1050, 2000, 5005, 1000],
        "id": ["c1", "c2", "c3", "c4", "c5"],
        "target_id": ["tc1", "tc2", "tc3", "tc4", "tc5"],
        "maf": [0.1, 0.5, 0.2, 0.3, 0.1],
        "priority": [2, 1, 1, 1, 1]
    })

def create_tsv(df, path):
    df.to_csv(path, sep="\t", index=False)

def test_basic_replacement(tmp_path, replace_data, candidate_data):
    replace_file = tmp_path / "replace.tsv"
    candidate_file = tmp_path / "candidate.tsv"
    out_file = tmp_path / "output.tsv"
    
    create_tsv(replace_data, replace_file)
    create_tsv(candidate_data, candidate_file)
    
    # Run basic command
    result = runner.invoke(app, [
        str(replace_file), 
        str(candidate_file), 
        str(out_file),
        "--window-kb", "1",
        "--max-distance-kb", "5"
    ])
    
    assert result.exit_code == 0
    assert out_file.exists()
    
    output_df = pd.read_table(out_file)
    # Check if 2 replacements were made
    assert len(output_df) == 2
    
    # Check logic for rs1 (pos 1000)
    # Candidates: 
    # c1: 1010 (dist 10), pri 2, maf 0.1
    # c2: 1050 (dist 50), pri 1, maf 0.5
    # c3: 2000 (dist 1000), pri 1, maf 0.2
    # Winner should be c2 because Priority 1 < 2.
    # Wait, logic says: Priority > Window > MAF > Distance > ID
    # c2 and c3 both Priority 1.
    # Window size 1kb (1000bp).
    # c2 distance 50 => window 0.
    # c3 distance 1000 => window 1.
    # So c2 wins on window index.
    
    row1 = output_df[output_df["origin_pos"] == 1000].iloc[0]
    assert row1["id"] == "c2"
    assert row1["origin_pos"] == 1000
    assert row1["distance_to_origin_bp"] == 50

    # Check logic for rs2 (pos 5000)
    # Candidates: c4 (5005, dist 5, pri 1).
    row2 = output_df[output_df["origin_pos"] == 5000].iloc[0]
    assert row2["id"] == "c4"


def test_no_candidates_strict(tmp_path, replace_data):
    replace_file = tmp_path / "replace.tsv"
    candidate_file = tmp_path / "candidate.tsv"
    out_file = tmp_path / "output.tsv"
    
    create_tsv(replace_data, replace_file)
    # Empty candidates on chr1
    candidate_df = pd.DataFrame({
        "chrom": ["2"],
        "pos": [1000],
        "id": ["c5"],
        "target_id": ["tc5"],
        "maf": [0.1],
        "priority": [1]
    })
    create_tsv(candidate_df, candidate_file)
    
    result = runner.invoke(app, [
        str(replace_file), 
        str(candidate_file), 
        str(out_file)
    ])
    
    assert result.exit_code != 0
    assert "找不到可用替换" in result.stdout or "错误" in result.stdout

def test_allow_missing(tmp_path):
    replace_file = tmp_path / "replace.tsv"
    candidate_file = tmp_path / "candidate.tsv"
    out_file = tmp_path / "output.tsv"
    
    # rs1 can be replaced, rs2 cannot
    replace_df = pd.DataFrame({
        "chrom": ["1", "1"],
        "pos": [1000, 90000],
        "id": ["rs1", "rs2"],
        "target_id": ["t1", "t2"]
    })
    
    candidate_df = pd.DataFrame({
        "chrom": ["1"],
        "pos": [1010],
        "id": ["c1"],
        "target_id": ["tc1"],
        "maf": [0.5],
        "priority": [1]
    })
    
    create_tsv(replace_df, replace_file)
    create_tsv(candidate_df, candidate_file)
    
    result = runner.invoke(app, [
        str(replace_file), 
        str(candidate_file), 
        str(out_file),
        "--allow-missing",
        "--max-distance-kb", "10" 
    ])
    
    assert result.exit_code == 0
    assert out_file.exists()
    
    # Check output contains only rs1 replacement
    output_df = pd.read_table(out_file)
    assert len(output_df) == 1
    assert output_df.iloc[0]["origin_pos"] == 1000
    
    # Check no_replacement file
    no_rep_file = tmp_path / "output_no_replacement.tsv"
    assert no_rep_file.exists()
    no_rep_df = pd.read_table(no_rep_file)
    assert len(no_rep_df) == 1
    assert no_rep_df.iloc[0]["id"] == "rs2"

def test_ranking_logic(tmp_path):
    replace_file = tmp_path / "replace.tsv"
    candidate_file = tmp_path / "candidate.tsv"
    out_file = tmp_path / "output.tsv"
    
    # 1 target
    replace_df = pd.DataFrame({
        "chrom": ["1"], "pos": [1000], "id": ["rs1"], "target_id": ["t1"]
    })
    
    # Candidates
    # c1: dist 10, pri 2, maf 0.9 (High MAF, close, but bad Priority)
    # c2: dist 20, pri 1, maf 0.1 (Good Priority, low MAF)
    # c3: dist 20, pri 1, maf 0.5 (Good Priority, higher MAF than c2)
    candidate_df = pd.DataFrame({
        "chrom": ["1", "1", "1"],
        "pos": [1010, 1020, 1020],
        "id": ["c1", "c2", "c3"],
        "target_id": ["t", "t", "t"],
        "maf": [0.9, 0.1, 0.5],
        "priority": [2, 1, 1]
    })
    
    create_tsv(replace_df, replace_file)
    create_tsv(candidate_df, candidate_file)
    
    result = runner.invoke(app, [str(replace_file), str(candidate_file), str(out_file)])
    
    assert result.exit_code == 0
    output_df = pd.read_table(out_file)
    # c3 wins: Priority 1 > 2 (beats c1). Between c2/c3 (both Pri 1, same dist/window), c3 has higher MAF.
    assert output_df.iloc[0]["id"] == "c3"


def test_ranking_window_logic(tmp_path):
    replace_file = tmp_path / "replace.tsv"
    candidate_file = tmp_path / "candidate.tsv"
    out_file = tmp_path / "output.tsv"
    
    # 1 target
    replace_df = pd.DataFrame({
        "chrom": ["1"], "pos": [1000], "id": ["rs1"], "target_id": ["t1"]
    })
    
    # Window 10kb.
    # c1: dist 100, pri 1, maf 0.1 (Window 0)
    # c2: dist 500, pri 1, maf 0.9 (Window 0) -> Winner due to MAF?
    # c3: dist 11000, pri 1, maf 0.99 (Window 1)
    
    # Wait, looking at _select_candidate:
    # sort_values(["priority", "window_index", "maf", "distance_bp", "id"], ascending=[True, True, False, True, True])
    
    # Within same window index, MAF is key.
    # c1 and c2 are in window 0 (0-10kb). c2 has higher MAF. c2 should win.
    # c3 is in window 1.
    
    candidate_df = pd.DataFrame({
        "chrom": ["1", "1", "1"],
        "pos": [1100, 1500, 12000],
        "id": ["c1", "c2", "c3"],
        "target_id": ["t", "t", "t"],
        "maf": [0.1, 0.9, 0.99],
        "priority": [1, 1, 1]
    })
    
    create_tsv(replace_df, replace_file)
    create_tsv(candidate_df, candidate_file)
    
    result = runner.invoke(app, [
        str(replace_file), 
        str(candidate_file), 
        str(out_file),
        "--window-kb", "10"
    ])
    
    assert result.exit_code == 0
    output_df = pd.read_table(out_file)
    assert output_df.iloc[0]["id"] == "c2"

def test_prevent_reuse(tmp_path):
    replace_file = tmp_path / "replace.tsv"
    candidate_file = tmp_path / "candidate.tsv"
    out_file = tmp_path / "output.tsv"
    
    # rs1 at 1000, rs2 at 1005. Very close.
    replace_df = pd.DataFrame({
        "chrom": ["1", "1"],
        "pos": [1000, 1005],
        "id": ["rs1", "rs2"],
        "target_id": ["t1", "t2"]
    })
    
    # Only one good candidate c1 at 1020.
    # c2 at 2000 (farther but valid).
    candidate_df = pd.DataFrame({
        "chrom": ["1", "1"],
        "pos": [1020, 2000],
        "id": ["c1", "c2"],
        "target_id": ["tc1", "tc2"],
        "maf": [0.5, 0.5],
        "priority": [1, 1]
    })
    
    create_tsv(replace_df, replace_file)
    create_tsv(candidate_df, candidate_file)
    
    result = runner.invoke(app, [str(replace_file), str(candidate_file), str(out_file)])
    assert result.exit_code == 0
    
    output_df = pd.read_table(out_file)
    ids = output_df["id"].tolist()
    # Should use c1 and c2. Shouldn't use c1 twice.
    assert "c1" in ids
    assert "c2" in ids
    assert len(set(ids)) == 2
