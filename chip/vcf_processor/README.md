# VCF Processor

An optimized, modular VCF processing system for genotyping arrays using cyvcf2.

## Two Usage Options

### Option 1: Standalone Script (Smart Hybrid)
```bash
# Process all variants in VCF file
python vcf_processor_standalone.py input.vcf output

# Process only specific target variants
python vcf_processor_standalone.py input.vcf output --targets targets.txt
```
- ✅ **Smart selection** - automatically uses high-performance modules if available, falls back to pure Python
- ✅ **Zero dependencies** - works with Python standard library only as fallback
- ✅ **No installation** - download and run immediately
- ✅ **Perfect for**: Quick processing, any environment, maximum compatibility

### Option 2: Full Module (High Performance)
```bash
# Process all variants in VCF file
python -m chip.vcf_processor.cli process input.vcf output

# Process only specific target variants
python -m chip.vcf_processor.cli process input.vcf output --targets targets.txt
```
- 🚀 **High performance** - uses cyvcf2 C extensions
- 🚀 **Full features** - configuration files, multi-threading, advanced options
- 🚀 **Perfect for**: Production environments, large-scale processing

**Both options use identical command-line arguments and produce identical outputs.**

## Overview

The VCF Processor is a high-performance tool designed to process VCF files for genotyping array workflows. It replaces bcftools-based processing with a pure Python implementation using cyvcf2 for improved performance and maintainability.

## Features

- **Fast VCF Processing**: Uses cyvcf2 for efficient VCF file parsing
- **Memory Efficient**: Configurable batch processing for large files
- **Flexible Output**: Supports both genotype (.gt.txt) and sequence (.seq.txt) formats
- **Variant Type Annotation**: Optional variant type classification (SNP/INDEL/MNP/REF)
- **Comprehensive CLI**: Full command-line interface with progress reporting
- **Error Handling**: Robust error handling with detailed logging
- **Performance Monitoring**: Built-in performance profiling and optimization
- **Extensive Testing**: Property-based testing with hypothesis for reliability

## Variant Type Annotation

The VCF Processor supports optional variant type annotation using the `--va-type` parameter. When enabled, it adds a `Variant_Type` column after the `ALT` column in both output files.

### Classification Rules

- **SNP**: Single nucleotide polymorphism (REF and ALT are both single bases, different)
- **INDEL**: Insertion or deletion (REF and ALT have different lengths)
- **MNP**: Multi-nucleotide polymorphism (REF and ALT have same length >1, different sequences)
- **REF**: Reference allele (REF and ALT are identical, no variation)

### Multi-Allelic Variants

For multi-allelic variants, types are combined with pipe separators and sorted in the order: `SNP|INDEL|MNP|REF`

Examples:
- Single SNP: `SNP`
- Single insertion: `INDEL`
- Multi-allelic with SNP and insertion: `SNP|INDEL`
- Multi-allelic with all types: `SNP|INDEL|MNP|REF`

### Usage Examples

```bash
# Enable variant type annotation
python vcf_processor_standalone.py input.vcf output --va-type

# With modular CLI
python -m chip.vcf_processor.cli process input.vcf output --va-type

# Combined with other options
python vcf_processor_standalone.py input.vcf output --targets targets.txt --va-type --compress
```

## Quick Start

### Option 1: Standalone Script (Recommended for Quick Start)

```bash
# Process all variants (no target file needed)
python vcf_processor_standalone.py input.vcf output

# Process specific target variants
python vcf_processor_standalone.py input.vcf output --targets targets.txt

# With variant type annotation
python vcf_processor_standalone.py input.vcf output --va-type

# With custom options (same as modular version)
python vcf_processor_standalone.py input.vcf output --targets targets.txt \
    --batch-size 20000 \
    --compress \
    --miss-fmt "NN" \
    --gt-sep "" \
    --va-type
```

**See [STANDALONE_USAGE.md](STANDALONE_USAGE.md) for detailed standalone usage guide.**

### Option 2: Full Module (Recommended for Production)

#### Installation

```bash
# Install dependencies
pip install cyvcf2 typer loguru tqdm pandas

# Or install from requirements
pip install -r requirements.txt
```

#### Usage

```bash
# Process all variants
python -m chip.vcf_processor.cli process input.vcf output

# Process specific target variants
python -m chip.vcf_processor.cli process input.vcf output --targets targets.txt

# With variant type annotation
python -m chip.vcf_processor.cli process input.vcf output --va-type

# With custom options (identical to standalone version)
python -m chip.vcf_processor.cli process input.vcf output --targets targets.txt \
    --batch-size 20000 \
    --threads 4 \
    --compress \
    --miss-fmt "NN" \
    --gt-sep "" \
    --va-type
```

### Python API

```python
from pathlib import Path
from chip.vcf_processor.config import ProcessingConfig
from chip.vcf_processor.vcf_processor import VCFProcessor

# Create configuration
config = ProcessingConfig(
    vcf_file=Path("input.vcf"),
    target_id_file=Path("targets.txt"),
    output_file=Path("output"),
    batch_size=10000,
    threads=4,
    compress_output=True,
    include_variant_type=True  # Enable variant type annotation
)

# Process VCF
processor = VCFProcessor(config)
result = processor.process()

# Get summary
summary = processor.format_summary()
print(summary)
```

## Choosing the Right Option

| Feature | Standalone Script | Full Module |
|---------|------------------|-------------|
| **Installation** | None required | `pip install` required |
| **Dependencies** | Auto-detects, falls back to stdlib | cyvcf2, typer, loguru, etc. |
| **Performance** | Good→Excellent (adaptive) | Excellent (C extensions) |
| **Memory Usage** | Efficient→More efficient | More efficient |
| **Multi-threading** | No | Yes (`--threads` option) |
| **Configuration Files** | No | Yes (JSON config support) |
| **Variant Type Annotation** | ✅ Yes (`--va-type`) | ✅ Yes (`--va-type`) |
| **Command Compatibility** | ✅ Identical arguments | ✅ Identical arguments |
| **Output Format** | ✅ Identical | ✅ Identical |
| **Best For** | Any environment, maximum compatibility | Production, large files |

**Recommendation**: The standalone script now intelligently adapts to your environment - it will use high-performance modules if available, or fall back to pure Python if needed. This gives you the best of both worlds!

## Architecture

The VCF Processor follows a modular architecture with separate components for each processing stage:

```
VCFProcessor (Orchestrator)
├── VCFReader (cyvcf2-based VCF reading)
├── VariantFilter (Target ID filtering)
├── VariantTransformer (Variant data transformation)
├── GenotypeConverter (Genotype format conversion)
├── OutputWriter (File output with compression)
└── ErrorHandler (Error management and logging)
```

### Key Components

- **VCFReader**: Efficient VCF file reading using cyvcf2
- **VariantFilter**: Filters variants based on target ID lists
- **VariantTransformer**: Transforms variant data for processing
- **GenotypeConverter**: Converts genotypes to desired output formats
- **OutputWriter**: Handles output file generation with optional compression
- **ProcessingConfig**: Centralized configuration management

## Module Structure

```
vcf_processor/
├── __init__.py                    # Package initialization
├── config.py                      # Configuration and data models
├── vcf_reader.py                 # VCF file reading with cyvcf2
├── variant_filter.py             # Target ID filtering
├── variant_transformer.py        # Variant data transformation
├── genotype_converter.py         # Genotype format conversion
├── output_writer.py              # File output management
├── vcf_processor.py              # Main orchestrator
├── performance_optimizer.py      # Performance optimization
├── error_handler.py              # Error management
├── logging_config.py             # Logging configuration
├── config_manager.py             # Configuration management
├── cli.py                        # Command-line interface
├── tests/                        # Comprehensive test suite
│   ├── test_*_unit.py           # Unit tests
│   ├── test_*_properties.py     # Property-based tests
│   ├── test_integration*.py     # Integration tests
│   └── conftest.py              # Test fixtures
└── README.md                     # This file
```

## Configuration

### Processing Configuration

```python
ProcessingConfig(
    vcf_file=Path("input.vcf"),           # Input VCF file
    target_id_file=Path("targets.txt"),   # Target variant IDs
    output_file=Path("output"),           # Output file prefix
    miss_fmt="NN",                        # Missing genotype format
    gt_sep="",                            # Genotype separator
    batch_size=10000,                     # Processing batch size
    threads=4,                            # Number of threads
    compress_output=True,                 # Enable output compression
    include_variant_type=False,           # Enable variant type annotation
    verbose=False,                        # Enable verbose logging
    quiet=False,                          # Enable quiet mode
    dry_run=False,                        # Preview mode
    log_file=None                         # Optional log file
)
```

### CLI Options

```bash
Options:
  --miss-fmt TEXT         Missing genotype format [default: NN]
  --gt-sep TEXT          Genotype separator [default: ]
  --threads INTEGER      Number of threads [default: 4]
  --batch-size INTEGER   Processing batch size [default: 10000]
  --compress/--no-compress  Compress output files [default: compress]
  --va-type              Add variant type annotation column
  --verbose / -v         Enable verbose logging
  --quiet / -q           Enable quiet mode
  --dry-run              Preview operations without executing
  --log-file PATH        Log file path
  --config PATH          Configuration file path
```

## Performance Optimization

### Batch Processing

For large files (>100MB), use batch processing:

```python
config = ProcessingConfig(
    batch_size=50000,  # Large batch size for memory efficiency
    threads=2,         # Moderate threading for I/O bound operations
    compress_output=True  # Save disk space
)
```

### Performance Tips

1. **Use appropriate batch sizes**:
   - Small files (<10MB): batch_size=1000-5000
   - Medium files (10-100MB): batch_size=10000-20000
   - Large files (>100MB): batch_size=50000+

2. **Enable compression for large outputs**:
   ```python
   config.compress_output = True
   ```

3. **Monitor memory usage**:
   ```python
   config.verbose = True  # Enable detailed logging
   ```

4. **Use dry-run for validation**:
   ```python
   config.dry_run = True  # Validate without processing
   ```

## Examples

### Example 1: Basic Processing

```python
from pathlib import Path
from chip.vcf_processor.config import ProcessingConfig
from chip.vcf_processor.vcf_processor import VCFProcessor

# Configure processing
config = ProcessingConfig(
    vcf_file=Path("data/genotypes.vcf"),
    target_id_file=Path("data/targets.txt"),
    output_file=Path("results/output")
)

# Process VCF
processor = VCFProcessor(config)
result = processor.process()

print(f"Successfully processed {result.processed_variants} variants")
```

### Example 2: Large File Processing

```python
# Optimized configuration for large files
config = ProcessingConfig(
    vcf_file=Path("data/large_genotypes.vcf.gz"),
    target_id_file=Path("data/targets.txt"),
    output_file=Path("results/large_output"),
    batch_size=50000,        # Large batches for efficiency
    threads=4,               # Parallel processing
    compress_output=True,    # Save disk space
    verbose=True            # Monitor progress
)

processor = VCFProcessor(config)
result = processor.process()

# Get detailed summary
summary = processor.format_summary()
print(summary)
```

### Example 3: Error Handling

```python
import logging
from chip.vcf_processor.config import ProcessingConfig
from chip.vcf_processor.vcf_processor import VCFProcessor

# Configure logging
logging.basicConfig(level=logging.INFO)

config = ProcessingConfig(
    vcf_file=Path("data/problematic.vcf"),
    target_id_file=Path("data/targets.txt"),
    output_file=Path("results/error_output"),
    verbose=True
)

try:
    processor = VCFProcessor(config)
    result = processor.process()
    
    if result.has_errors:
        print("Processing completed with errors:")
        for error in result.errors:
            print(f"  - {error}")
    else:
        print("Processing completed successfully")
        
except Exception as e:
    print(f"Processing failed: {e}")
```

## Testing

The VCF Processor includes comprehensive testing:

### Running Tests

```bash
# Run all tests
pytest chip/vcf_processor/tests/

# Run with coverage
pytest --cov=chip.vcf_processor chip/vcf_processor/tests/

# Run property-based tests
pytest chip/vcf_processor/tests/test_*_properties.py

# Run integration tests
pytest chip/vcf_processor/tests/test_integration*.py

# Run specific test categories
pytest chip/vcf_processor/tests/ -k "unit"
pytest chip/vcf_processor/tests/ -k "property"
pytest chip/vcf_processor/tests/ -k "integration"
```

### Test Categories

- **Unit Tests**: Test individual components in isolation
- **Property Tests**: Test universal properties with hypothesis
- **Integration Tests**: Test complete workflows end-to-end
- **Performance Tests**: Benchmark processing performance

### Test Coverage

The test suite maintains >90% code coverage and includes:

- 17 property-based tests validating universal correctness properties
- Comprehensive unit tests for all components
- Integration tests with various VCF formats and sizes
- Error scenario testing and recovery validation
- Performance benchmarking and optimization validation

## Migration Guide

### From bcftools-based Processing

1. **Update dependencies**:
   ```bash
   pip install cyvcf2 typer loguru tqdm
   ```

2. **Update command-line usage**:
   ```bash
   # Old bcftools approach
   bcftools view input.vcf | process_script.py
   
   # New VCF processor
   python -m chip.vcf_processor.cli process input.vcf targets.txt output
   ```

3. **Update Python code**:
   ```python
   # Old approach
   import subprocess
   result = subprocess.run(["bcftools", "view", "input.vcf"])
   
   # New approach
   from chip.vcf_processor.vcf_processor import VCFProcessor
   processor = VCFProcessor(config)
   result = processor.process()
   ```

## Troubleshooting

### Common Issues

1. **ImportError: No module named 'cyvcf2'**
   ```bash
   pip install cyvcf2
   ```

2. **Memory errors with large files**
   - Reduce batch_size: `config.batch_size = 5000`
   - Enable compression: `config.compress_output = True`
   - Close other applications to free memory

3. **Slow processing**
   - Increase batch_size for large files: `config.batch_size = 50000`
   - Use appropriate thread count: `config.threads = 4`
   - Enable progress monitoring: `config.verbose = True`

4. **Permission denied errors**
   - Check output directory permissions
   - Ensure sufficient disk space
   - Use different output directory

## API Reference

### Main Classes

- **ProcessingConfig**: Configuration management
- **VCFProcessor**: Main processing orchestrator
- **VCFProcessorFactory**: Factory for creating processors
- **ProcessingResult**: Container for processing results

### Key Methods

- `VCFProcessor.process()`: Execute VCF processing
- `VCFProcessor.validate_output()`: Validate output files
- `VCFProcessor.get_processing_summary()`: Get processing statistics
- `VCFProcessor.format_summary()`: Get formatted summary

## Contributing

### Development Setup

```bash
# Clone repository
git clone <repository-url>
cd tc-pytools

# Install development dependencies
pip install -r requirements.txt
pip install -e .

# Install pre-commit hooks
pre-commit install

# Run tests
pytest chip/vcf_processor/tests/
```

### Code Style

- Follow Google Python Style Guide
- Use type hints for all functions
- Add docstrings for all public methods
- Format code with ruff: `ruff format .`
- Lint code with ruff: `ruff check . --fix`

### Testing Requirements

- Add unit tests for new functionality
- Include property-based tests for complex logic
- Add integration tests for end-to-end workflows
- Ensure >90% test coverage

## License

This project is part of TC PyTools and follows the same licensing terms.