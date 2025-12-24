# TC PyTools Product Overview

TC PyTools is a collection of bioinformatics Python scripts and tools organized by business scenarios (panel, primer, fasta, vcf, etc.). This repository serves as a toolkit rather than a unified application - each script operates independently with its own purpose and documentation.

## Core Purpose

- **Bioinformatics Pipeline Tools**: Scripts for common genomics workflows including VCF processing, primer design, FASTA manipulation, and panel analysis
- **Modular Design**: Each module (chip, panel, primer, etc.) contains related functionality with independent scripts
- **Research & Production**: Tools used in both research and production genomics pipelines

## Key Modules

- **chip**: Genotyping array and VCF processing tools
- **panel**: Target panel design and coverage analysis
- **primer**: Primer design and KASP assay tools  
- **vcf**: VCF file processing and statistics
- **fasta**: FASTA sequence manipulation
- **gwas**: GWAS analysis and visualization
- **exome**: Exome sequencing analysis tools

## Usage Pattern

Scripts are designed to be run independently from their respective directories. Each module contains its own README with specific usage instructions and examples.