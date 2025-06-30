#!/usr/bin/env python3
"""
Test script for add_flank_sequence.py
"""

import os
import tempfile

import pandas as pd

# Create a sample table
sample_data = {
    "chrom": ["chr1", "chr1", "chr2"],
    "pos": [1000, 2000, 1500],  # 1-based positions
    "probe_start": [995, 1995, 1495],  # 0-based
    "probe_end": [1005, 2005, 1505],  # 0-based
    "alleles": ["A/T", "G/C", "T/A"],
}

df = pd.DataFrame(sample_data)

# Save sample table
with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
    df.to_csv(f.name, sep="\t", index=False)
    table_file = f.name

print(f"Sample table created: {table_file}")
print("Table contents:")
print(df.to_string())

print("\nTo use the script:")
print(
    f"python add_flank_sequence.py {table_file} <genome.fa> output.tsv --flank-size 20"
)

# Clean up
os.unlink(table_file)
