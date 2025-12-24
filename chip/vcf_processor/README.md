# VCF Processor Optimization

This module provides an optimized VCF genotype table processor that replaces bcftools with cyvcf2 for better performance and maintainability.

## Features

- **cyvcf2 Integration**: Fast VCF parsing without external dependencies
- **Modular Architecture**: Clean separation of concerns with dedicated classes
- **Comprehensive Testing**: Unit tests and property-based tests
- **Error Handling**: Robust error recovery and detailed error reporting
- **Performance Optimized**: Streaming processing and configurable batch sizes
- **Backward Compatible**: Maintains identical output format to original implementation

## Module Structure

```
vcf_processor/
├── __init__.py              # Package initialization
├── config.py                # Configuration and data models
├── vcf_reader.py           # VCF file reading with cyvcf2
├── variant_filter.py       # Target ID filtering
├── genotype_converter.py   # Genotype to sequence conversion
├── variant_transformer.py  # Variant notation transformation
├── output_writer.py        # File output management
├── processor.py            # Main orchestrator
├── logging_config.py       # Logging setup
├── cli.py                  # Command-line interface
├── tests/                  # Test suite
└── README.md              # This file
```

## Usage

```python
from chip.vcf_processor import VCFProcessor, ProcessingConfig

# Create configuration
config = ProcessingConfig(
    vcf_file=Path("input.vcf.gz"),
    target_id_file=Path("targets.txt"),
    output_file=Path("output"),
    miss_fmt="NN",
    threads=4
)

# Process VCF
processor = VCFProcessor(config)
result = processor.process()
print(result.summary())
```

## Command Line

```bash
python -m chip.vcf_processor.cli input.vcf.gz targets.txt output --threads 4
```

## Testing

```bash
# Run all tests
pytest chip/vcf_processor/tests/

# Run with coverage
pytest chip/vcf_processor/tests/ --cov=chip.vcf_processor

# Run property-based tests only
pytest chip/vcf_processor/tests/ -k "property"
```