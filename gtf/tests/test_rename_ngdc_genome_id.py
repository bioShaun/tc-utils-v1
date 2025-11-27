import tempfile
import unittest
from pathlib import Path

from gtf.rename_ngdc_genome_id import (
    build_id_mapping,
    parse_fasta_header,
    rename_fasta,
    rename_gff,
)


class TestParseFastaHeader(unittest.TestCase):
    """Test parsing of FASTA header to extract old and new IDs."""
    
    def test_standard_ngdc_format(self):
        header = ">GWHGECT00000001.1      Chromosome 1A   Complete=T      Circular=F      OriSeqID=Chr1A  Len=600907804"
        old_id, new_id = parse_fasta_header(header)
        self.assertEqual(old_id, "GWHGECT00000001.1")
        self.assertEqual(new_id, "Chr1A")
    
    def test_another_example(self):
        header = ">GWHGECT00000010.1      Chromosome 5B   Complete=T      Circular=F      OriSeqID=Chr5B  Len=742642744"
        old_id, new_id = parse_fasta_header(header)
        self.assertEqual(old_id, "GWHGECT00000010.1")
        self.assertEqual(new_id, "Chr5B")
    
    def test_header_without_oriseqid(self):
        header = ">GWHGECT00000001.1      Chromosome 1A"
        old_id, new_id = parse_fasta_header(header)
        self.assertIsNone(old_id)
        self.assertIsNone(new_id)
    
    def test_empty_header(self):
        header = ">"
        old_id, new_id = parse_fasta_header(header)
        self.assertIsNone(old_id)
        self.assertIsNone(new_id)


class TestBuildIdMapping(unittest.TestCase):
    """Test building ID mapping from FASTA file."""
    
    def test_build_mapping(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            fasta_content = """>GWHGECT00000001.1      Chromosome 1A   Complete=T      Circular=F      OriSeqID=Chr1A  Len=600907804
ATCGATCGATCG
>GWHGECT00000002.1      Chromosome 1B   Complete=T      Circular=F      OriSeqID=Chr1B  Len=731628012
GCTAGCTAGCTA
>GWHGECT00000003.1      Chromosome 2A   Complete=T      Circular=F      OriSeqID=Chr2A  Len=801619444
TTAATTAATTAA
"""
            fasta_file = Path(tmp_dir) / "test.fasta"
            fasta_file.write_text(fasta_content)
            
            id_map = build_id_mapping(fasta_file)
            
            self.assertEqual(len(id_map), 3)
            self.assertEqual(id_map["GWHGECT00000001.1"], "Chr1A")
            self.assertEqual(id_map["GWHGECT00000002.1"], "Chr1B")
            self.assertEqual(id_map["GWHGECT00000003.1"], "Chr2A")


class TestRenameFasta(unittest.TestCase):
    """Test renaming chromosome IDs in FASTA file."""
    
    def test_rename_fasta(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            input_content = """>GWHGECT00000001.1      Chromosome 1A   Complete=T      Circular=F      OriSeqID=Chr1A  Len=600907804
ATCGATCGATCG
GCTAGCTAGCTA
>GWHGECT00000002.1      Chromosome 1B   Complete=T      Circular=F      OriSeqID=Chr1B  Len=731628012
TTAATTAATTAA
CCGGCCGGCCGG
"""
            input_file = Path(tmp_dir) / "input.fasta"
            output_file = Path(tmp_dir) / "output.fasta"
            input_file.write_text(input_content)
            
            id_map = {
                "GWHGECT00000001.1": "Chr1A",
                "GWHGECT00000002.1": "Chr1B"
            }
            
            rename_fasta(input_file, output_file, id_map)
            
            output_content = output_file.read_text()
            expected = """>Chr1A
ATCGATCGATCG
GCTAGCTAGCTA
>Chr1B
TTAATTAATTAA
CCGGCCGGCCGG
"""
            self.assertEqual(output_content, expected)
    
    def test_rename_fasta_no_mapping(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            input_content = """>GWHGECT00000001.1      Chromosome 1A   Complete=T      Circular=F      OriSeqID=Chr1A  Len=600907804
ATCG
>UnknownSeq
GCTA
"""
            input_file = Path(tmp_dir) / "input.fasta"
            output_file = Path(tmp_dir) / "output.fasta"
            input_file.write_text(input_content)
            
            id_map = {"GWHGECT00000001.1": "Chr1A"}
            
            rename_fasta(input_file, output_file, id_map)
            
            output_content = output_file.read_text()
            expected = """>Chr1A
ATCG
>UnknownSeq
GCTA
"""
            self.assertEqual(output_content, expected)


class TestRenameGff(unittest.TestCase):
    """Test renaming chromosome IDs in GFF file."""
    
    def test_rename_gff(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            input_content = """##gff-version 3
##sequence-region GWHGECT00000001.1 1 600907804
GWHGECT00000001.1\tRefSeq\tgene\t100\t200\t.\t+\t.\tID=gene1;Name=GENE1
GWHGECT00000001.1\tRefSeq\tmRNA\t100\t200\t.\t+\t.\tID=transcript1;Parent=gene1
GWHGECT00000002.1\tRefSeq\tgene\t300\t400\t.\t-\t.\tID=gene2;Name=GENE2
"""
            input_file = Path(tmp_dir) / "input.gff"
            output_file = Path(tmp_dir) / "output.gff"
            input_file.write_text(input_content)
            
            id_map = {
                "GWHGECT00000001.1": "Chr1A",
                "GWHGECT00000002.1": "Chr1B"
            }
            
            rename_gff(input_file, output_file, id_map)
            
            output_content = output_file.read_text()
            expected = """##gff-version 3
##sequence-region GWHGECT00000001.1 1 600907804
Chr1A\tRefSeq\tgene\t100\t200\t.\t+\t.\tID=gene1;Name=GENE1
Chr1A\tRefSeq\tmRNA\t100\t200\t.\t+\t.\tID=transcript1;Parent=gene1
Chr1B\tRefSeq\tgene\t300\t400\t.\t-\t.\tID=gene2;Name=GENE2
"""
            self.assertEqual(output_content, expected)
    
    def test_rename_gff_with_unmapped_chromosomes(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            input_content = """GWHGECT00000001.1\tRefSeq\tgene\t100\t200\t.\t+\t.\tID=gene1
UnknownChr\tRefSeq\tgene\t300\t400\t.\t-\t.\tID=gene2
GWHGECT00000002.1\tRefSeq\tgene\t500\t600\t.\t+\t.\tID=gene3
"""
            input_file = Path(tmp_dir) / "input.gff"
            output_file = Path(tmp_dir) / "output.gff"
            input_file.write_text(input_content)
            
            id_map = {
                "GWHGECT00000001.1": "Chr1A",
                "GWHGECT00000002.1": "Chr1B"
            }
            
            rename_gff(input_file, output_file, id_map)
            
            output_content = output_file.read_text()
            expected = """Chr1A\tRefSeq\tgene\t100\t200\t.\t+\t.\tID=gene1
UnknownChr\tRefSeq\tgene\t300\t400\t.\t-\t.\tID=gene2
Chr1B\tRefSeq\tgene\t500\t600\t.\t+\t.\tID=gene3
"""
            self.assertEqual(output_content, expected)
    
    def test_rename_gff_empty_lines(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            input_content = """##gff-version 3

GWHGECT00000001.1\tRefSeq\tgene\t100\t200\t.\t+\t.\tID=gene1

# This is a comment
GWHGECT00000002.1\tRefSeq\tgene\t300\t400\t.\t-\t.\tID=gene2
"""
            input_file = Path(tmp_dir) / "input.gff"
            output_file = Path(tmp_dir) / "output.gff"
            input_file.write_text(input_content)
            
            id_map = {
                "GWHGECT00000001.1": "Chr1A",
                "GWHGECT00000002.1": "Chr1B"
            }
            
            rename_gff(input_file, output_file, id_map)
            
            output_content = output_file.read_text()
            expected = """##gff-version 3

Chr1A\tRefSeq\tgene\t100\t200\t.\t+\t.\tID=gene1

# This is a comment
Chr1B\tRefSeq\tgene\t300\t400\t.\t-\t.\tID=gene2
"""
            self.assertEqual(output_content, expected)


if __name__ == '__main__':
    unittest.main()
