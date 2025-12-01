import pytest
import pandas as pd
from pathlib import Path
from utils.tsv2excel import file2excel

def test_tsv_to_excel(tmp_path):
    # Create a dummy TSV file
    tsv_file = tmp_path / "test.tsv"
    tsv_file.write_text("col1\tcol2\n1\t2\n3\t4", encoding="utf-8")
    
    excel_file = tmp_path / "test.xlsx"
    
    shape = file2excel(tsv_file, excel_file)
    
    assert shape == (2, 2)
    assert excel_file.exists()
    
    df = pd.read_excel(excel_file)
    assert df.shape == (2, 2)
    assert df.iloc[0, 0] == 1

def test_csv_to_excel(tmp_path):
    # Create a dummy CSV file
    csv_file = tmp_path / "test.csv"
    csv_file.write_text("col1,col2\n1,2\n3,4", encoding="utf-8")
    
    excel_file = tmp_path / "test.xlsx"
    
    shape = file2excel(csv_file, excel_file)
    
    assert shape == (2, 2)
    assert excel_file.exists()
    
    df = pd.read_excel(excel_file)
    assert df.shape == (2, 2)

def test_auto_detect_tsv(tmp_path):
    # Create a dummy file without extension but tab separated
    txt_file = tmp_path / "test.txt"
    txt_file.write_text("col1\tcol2\n1\t2", encoding="utf-8")
    
    excel_file = tmp_path / "test.xlsx"
    
    shape = file2excel(txt_file, excel_file)
    assert shape == (1, 2)

def test_file_not_found():
    with pytest.raises(FileNotFoundError):
        file2excel(Path("non_existent.tsv"), Path("output.xlsx"))
