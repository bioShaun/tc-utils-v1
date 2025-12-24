# Technology Stack & Build System

## Core Technologies

- **Python 3.11+**: Primary language for all scripts
- **Package Management**: Uses `pyproject.toml` with setuptools backend
- **Dependency Management**: `requirements.txt` and `pyproject.toml` dependencies

## Key Libraries

### Data Processing
- **pandas**: Data manipulation and analysis
- **polars**: High-performance DataFrame library for large datasets
- **numpy**: Numerical computing

### Bioinformatics
- **pysam**: SAM/BAM file processing
- **cyvcf2**: Fast VCF file parsing (preferred over bcftools)
- **pyfaidx**: FASTA file indexing and access

### CLI & Utilities
- **typer**: Command-line interface framework
- **tqdm**: Progress bars
- **loguru**: Structured logging
- **delegator.py**: Process execution (with subprocess fallback)

### Validation & Testing
- **pandera**: Data validation
- **attrs**: Class definitions
- **pytest**: Testing framework
- **hypothesis**: Property-based testing

## Development Tools

### Code Quality
- **ruff**: Linting and formatting (replaces black, isort, flake8)
- **pyright**: Type checking
- **pre-commit**: Git hooks for automated checks

### Testing
- **pytest**: Primary testing framework
- **pytest-cov**: Coverage reporting
- **hypothesis**: Property-based testing for robust validation

## Common Commands

### Setup & Installation
```bash
# Install dependencies
pip install -r requirements.txt

# Install in editable mode (for development)
pip install -e .

# Setup pre-commit hooks
pre-commit install
```

### Development Workflow
```bash
# Run all quality checks
pre-commit run --all-files

# Run tests
pytest

# Run tests with coverage
pytest --cov=. --cov-report=term-missing

# Run specific module tests
pytest panel/tests/ -v
```

### Code Quality
```bash
# Format code
ruff format .

# Lint and auto-fix
ruff check . --fix

# Type checking
pyright
```

## Build Configuration

- **Line Length**: 120 characters
- **Python Version**: 3.11+ required
- **Import Organization**: Handled by ruff (isort replacement)
- **Type Checking**: Basic mode with pyright
- **Test Discovery**: Automatic via pytest configuration