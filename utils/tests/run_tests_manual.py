import unittest
import pandas as pd
from pathlib import Path
import shutil
import tempfile
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from utils.tsv2excel import file2excel

class TestTsv2Excel(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.test_path = Path(self.test_dir)

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_tsv_to_excel(self):
        tsv_file = self.test_path / "test.tsv"
        with open(tsv_file, "w", encoding="utf-8") as f:
            f.write("col1\tcol2\n1\t2\n3\t4")
        
        excel_file = self.test_path / "test.xlsx"
        
        shape = file2excel(tsv_file, excel_file)
        
        self.assertEqual(shape, (2, 2))
        self.assertTrue(excel_file.exists())
        
        df = pd.read_excel(excel_file)
        self.assertEqual(df.shape, (2, 2))
        self.assertEqual(df.iloc[0, 0], 1)

    def test_csv_to_excel(self):
        csv_file = self.test_path / "test.csv"
        with open(csv_file, "w", encoding="utf-8") as f:
            f.write("col1,col2\n1,2\n3,4")
        
        excel_file = self.test_path / "test.xlsx"
        
        shape = file2excel(csv_file, excel_file)
        
        self.assertEqual(shape, (2, 2))
        self.assertTrue(excel_file.exists())
        
        df = pd.read_excel(excel_file)
        self.assertEqual(df.shape, (2, 2))

    def test_auto_detect_tsv(self):
        txt_file = self.test_path / "test.txt"
        with open(txt_file, "w", encoding="utf-8") as f:
            f.write("col1\tcol2\n1\t2")
        
        excel_file = self.test_path / "test.xlsx"
        
        shape = file2excel(txt_file, excel_file)
        self.assertEqual(shape, (1, 2))

    def test_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            file2excel(Path("non_existent.tsv"), Path("output.xlsx"))

if __name__ == "__main__":
    unittest.main()
