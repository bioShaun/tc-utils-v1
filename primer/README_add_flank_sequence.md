# Add Probe Sequence Script

This script adds probe sequences to a table containing genomic positions and alleles.

## Usage

```bash
python add_flank_sequence.py <table_file> <genome_fasta> <output_file>
```

## Parameters

- `table_file`: Input table file (CSV or TSV format)
- `genome_fasta`: Genome FASTA file
- `output_file`: Output table file

## Input Table Format

The input table must contain the following columns:

- `chrom`: Chromosome name (e.g., "chr1", "chr2")
- `pos`: 1-based genomic position
- `probe_start`: Probe start position (0-based)
- `probe_end`: Probe end position (0-based)
- `alleles`: Alleles in format like "A/T", "G/C", "T/A"

## Example Input Table

```tsv
chrom	pos	probe_start	probe_end	alleles
chr1	1000	995	1005	A/T
chr1	2000	1995	2005	G/C
chr2	1500	1495	1505	T/A
```

## Output

The script adds a new `probe_sequence` column containing sequences where:
- The sequence is extracted from the genome from `probe_start` to `probe_end`
- The position at `pos` is replaced with the alleles in `[A/T]` format
- The sequence length is determined by `probe_end - probe_start`

## Example Output

```tsv
chrom	pos	probe_start	probe_end	alleles	probe_sequence
chr1	1000	995	1005	A/T	GCT[A/T]GAT
chr1	2000	1995	2005	G/C	TAG[G/C]CTA
chr2	1500	1495	1505	T/A	ATC[T/A]GCT
```

## How it works

1. For each row, extract sequence from `probe_start` to `probe_end` (0-based coordinates)
2. Calculate the position within the probe sequence: `pos_in_probe = pos - 1 - probe_start`
3. Replace the base at that position with `[alleles]` format
4. Add the modified sequence to the `probe_sequence` column

## Dependencies

- pandas
- pysam
- typer

## Example Usage

```bash
# Basic usage
python add_flank_sequence.py input.tsv genome.fa output.tsv

# Using CSV input
python add_flank_sequence.py input.csv genome.fa output.csv
``` 