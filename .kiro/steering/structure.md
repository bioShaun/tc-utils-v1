# Project Structure & Organization

## Repository Layout

TC PyTools follows a **modular structure** where each top-level directory represents a functional domain in bioinformatics workflows.

### Core Module Structure
```
module_name/
├── __init__.py              # Package initialization
├── script1.py               # Independent executable scripts
├── script2.py
├── README.md                # Module-specific documentation
└── tests/                   # Module test suite
    ├── __init__.py
    ├── test_script1.py
    └── conftest.py          # pytest fixtures
```

### Key Modules

- **`bsa/`**: Bulk Segregant Analysis tools
- **`chip/`**: Genotyping array and SNP chip processing
  - `vcf_processor/`: Modular VCF processing system (example of structured submodule)
- **`database/`**: Database interaction scripts
- **`exome/`**: Exome sequencing analysis tools
- **`fasta/`**: FASTA sequence manipulation utilities
- **`fq/`**: FASTQ file processing tools
- **`gtf/`**: GTF/GFF annotation processing
- **`gwas/`**: GWAS analysis and visualization
- **`liftover/`**: Genome coordinate conversion
- **`panel/`**: Target panel design and analysis
- **`primer/`**: Primer design and validation
- **`vcf/`**: VCF file processing utilities
- **`utils/`**: General utility functions

### Configuration & Metadata

- **`.kiro/`**: Kiro IDE configuration and steering rules
- **`conductor/`**: Project management and guidelines
- **`test-data/`**: Shared test datasets
- **`*-notes/`**: AI-generated documentation (codex, gemini, minimax)

## Architectural Patterns

### Script Independence
- Each script is **self-contained** and executable independently
- No shared entry point - scripts run from their module directories
- Module-level `__init__.py` files for package structure but minimal interdependence

### Modular Design Principles
- **Domain Separation**: Each module focuses on specific bioinformatics workflows
- **Consistent Structure**: All modules follow similar organization patterns
- **Documentation Co-location**: README files alongside code in each module

### Advanced Modules
Some modules like `chip/vcf_processor/` demonstrate **structured architecture**:
- Separate classes for different responsibilities (reader, filter, transformer, writer)
- Configuration management with dataclasses
- Comprehensive test suites with property-based testing
- CLI interfaces using typer

### Testing Strategy
- **Module-level tests**: Each module has its own `tests/` directory
- **pytest configuration**: Centralized in root `pytest.ini` with module-specific test paths
- **Property-based testing**: Used for robust validation (especially in vcf_processor)
- **Coverage reporting**: Integrated with pytest-cov

### File Naming Conventions
- **Scripts**: Descriptive names with hyphens (e.g., `align-kasp.py`, `simple-vcf-stats.py`)
- **Modules**: Snake_case for Python modules (e.g., `vcf_processor`, `gene_expression`)
- **Tests**: Prefix with `test_` following pytest conventions
- **Documentation**: `README.md` in each module, specialized docs in `*-notes/` directories

### Dependencies & Imports
- **External dependencies**: Managed via `requirements.txt` and `pyproject.toml`
- **Internal imports**: Minimal cross-module dependencies
- **Optional dependencies**: Graceful fallbacks (e.g., delegator.py → subprocess)