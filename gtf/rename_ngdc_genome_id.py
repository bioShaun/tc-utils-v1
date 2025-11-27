#!/usr/bin/env python3
"""
Rename chromosome IDs in NGDC genome FASTA and GFF files.
Extracts OriSeqID from FASTA headers and replaces GWHGECT IDs with the OriSeqID names.

Usage:
    python rename_ngdc_genome_id.py -f genome.fasta -o output.fasta [-g input.gff -og output.gff]
"""

import argparse
import re
import sys


def parse_fasta_header(header):
    """
    Parse FASTA header to extract ID and OriSeqID.
    
    Example header:
    >GWHGECT00000001.1      Chromosome 1A   Complete=T      Circular=F      OriSeqID=Chr1A  Len=600907804
    
    Returns:
        tuple: (old_id, new_id) or (None, None) if OriSeqID not found
    """
    # Extract the first field (ID)
    parts = header.strip().split()
    if not parts:
        return None, None
    
    old_id = parts[0].lstrip('>')
    
    # Extract OriSeqID
    match = re.search(r'OriSeqID=(\S+)', header)
    if match:
        new_id = match.group(1)
        return old_id, new_id
    
    return None, None


def build_id_mapping(fasta_file):
    """
    Build a mapping dictionary from FASTA file.
    
    Args:
        fasta_file: Path to input FASTA file
        
    Returns:
        dict: Mapping from old IDs to new IDs
    """
    id_map = {}
    
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                old_id, new_id = parse_fasta_header(line)
                if old_id and new_id:
                    id_map[old_id] = new_id
    
    return id_map


def rename_fasta(input_fasta, output_fasta, id_map):
    """
    Rename chromosome IDs in FASTA file.
    
    Args:
        input_fasta: Path to input FASTA file
        output_fasta: Path to output FASTA file
        id_map: Dictionary mapping old IDs to new IDs
    """
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                old_id, new_id = parse_fasta_header(line)
                if old_id and new_id and old_id in id_map:
                    # Write simplified header with just the new ID
                    outfile.write(f'>{new_id}\n')
                else:
                    # Keep original header if no mapping found
                    outfile.write(line)
            else:
                outfile.write(line)


def rename_gff(input_gff, output_gff, id_map):
    """
    Rename chromosome IDs in GFF file.
    
    Args:
        input_gff: Path to input GFF file
        output_gff: Path to output GFF file
        id_map: Dictionary mapping old IDs to new IDs
    """
    with open(input_gff, 'r') as infile, open(output_gff, 'w') as outfile:
        for line in infile:
            # Skip comment lines
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            # Skip empty lines
            if not line.strip():
                outfile.write(line)
                continue
            
            # Process GFF lines
            fields = line.split('\t')
            if len(fields) >= 1:
                chrom = fields[0]
                if chrom in id_map:
                    fields[0] = id_map[chrom]
                outfile.write('\t'.join(fields))
            else:
                outfile.write(line)


def main():
    parser = argparse.ArgumentParser(
        description='Rename chromosome IDs in NGDC genome FASTA and GFF files using OriSeqID'
    )
    parser.add_argument('-f', '--fasta', required=True,
                        help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output FASTA file')
    parser.add_argument('-g', '--gff', required=False,
                        help='Input GFF file (optional)')
    parser.add_argument('-og', '--output-gff', required=False,
                        help='Output GFF file (optional, required if -g is specified)')
    
    args = parser.parse_args()
    
    # Check GFF arguments consistency
    if args.gff and not args.output_gff:
        parser.error('--output-gff is required when --gff is specified')
    if args.output_gff and not args.gff:
        parser.error('--gff is required when --output-gff is specified')
    
    # Build ID mapping from FASTA
    print(f'Building ID mapping from {args.fasta}...', file=sys.stderr)
    id_map = build_id_mapping(args.fasta)
    print(f'Found {len(id_map)} chromosome mappings', file=sys.stderr)
    
    for old_id, new_id in sorted(id_map.items()):
        print(f'  {old_id} -> {new_id}', file=sys.stderr)
    
    # Rename FASTA
    print(f'Renaming FASTA file to {args.output}...', file=sys.stderr)
    rename_fasta(args.fasta, args.output, id_map)
    
    # Rename GFF if specified
    if args.gff:
        print(f'Renaming GFF file to {args.output_gff}...', file=sys.stderr)
        rename_gff(args.gff, args.output_gff, id_map)
    
    print('Done!', file=sys.stderr)


if __name__ == '__main__':
    main()
