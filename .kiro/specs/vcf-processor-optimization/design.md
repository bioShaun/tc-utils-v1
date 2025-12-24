# Design Document: VCF Processor Optimization

## Overview

This design outlines the optimization of the existing VCF genotype table processing script by replacing bcftools with cyvcf2, improving code structure, and adding comprehensive testing. The new implementation will maintain backward compatibility while providing better performance, maintainability, and reliability.

## Architecture

The optimized VCF processor will follow a modular architecture with clear separation of concerns:

```
VCFProcessor (Main Orchestrator)
├── VCFReader (cyvcf2-based VCF parsing)
├── VariantFilter (Target ID filtering)
├── GenotypeConverter (GT to sequence conversion)
├── VariantTransformer (VCF notation to annotation format)
├── OutputWriter (File output management)
└── ConfigManager (Configuration and CLI handling)
```

## Components and Interfaces

### VCFReader Class
```python
class VCFReader:
    def __init__(self, vcf_path: Path)
    def get_sample_names(self) -> List[str]
    def iter_variants(self, target_ids: Optional[Set[str]] = None) -> Iterator[Variant]
    def get_variant_by_id(self, variant_id: str) -> Optional[Variant]
```

### VariantFilter Class
```python
class VariantFilter:
    def __init__(self, target_ids: Set[str])
    def should_include(self, variant: Variant) -> bool
    def load_target_ids(self, target_file: Path) -> Set[str]
```

### GenotypeConverter Class
```python
class GenotypeConverter:
    def __init__(self, miss_fmt: str = "NN", gt_sep: str = "")
    def convert_genotype(self, variant: Variant, sample_idx: int) -> str
    def convert_batch(self, variants: List[Variant]) -> pd.DataFrame
```

### VariantTransformer Class
```python
class VariantTransformer:
    @staticmethod
    def transform_alt(ref: str, alt: str) -> str
    @staticmethod
    def transform_multi_alt(ref: str, alts: List[str]) -> str
    def validate_transformation(self, ref: str, alt: str) -> bool
```

### OutputWriter Class
```python
class OutputWriter:
    def __init__(self, output_path: Path, compress: bool = True)
    def write_genotype_table(self, df: pd.DataFrame, append: bool = False)
    def write_sequence_table(self, df: pd.DataFrame, append: bool = False)
    def __enter__(self) / __exit__(self) # Context manager support
```

### ConfigManager Class
```python
@dataclass
class ProcessingConfig:
    vcf_file: Path
    target_id_file: Path
    output_file: Path
    miss_fmt: str = "NN"
    gt_sep: str = ""
    threads: int = 4
    batch_size: int = 10000
    compress_output: bool = True
    verbose: bool = False

class ConfigManager:
    @staticmethod
    def from_cli() -> ProcessingConfig
    @staticmethod
    def from_file(config_path: Path) -> ProcessingConfig
    def validate(self) -> None
```

## Data Models

### Variant Data Structure
```python
@dataclass
class VariantInfo:
    chrom: str
    pos: int
    ref: str
    alt: List[str]
    variant_id: str
    genotypes: List[str]  # Raw GT values
    
    def to_dict(self) -> Dict[str, Any]
    def get_location_key(self) -> str
```

### Processing Result
```python
@dataclass
class ProcessingResult:
    total_variants: int
    processed_variants: int
    skipped_variants: int
    errors: List[str]
    processing_time: float
    
    def summary(self) -> str
```

## Error Handling

The system will implement comprehensive error handling:

1. **VCF File Errors**: Invalid format, corrupted files, missing samples
2. **Target ID Errors**: Invalid format, non-existent IDs
3. **Transformation Errors**: Complex variants that cannot be transformed
4. **I/O Errors**: File permission issues, disk space problems
5. **Memory Errors**: Large file handling with graceful degradation

Each error type will have specific handling strategies:
- Log detailed error information
- Continue processing when possible
- Provide recovery suggestions
- Generate error reports

## Testing Strategy

### Unit Testing
- Test each component in isolation
- Mock external dependencies (file I/O, cyvcf2 objects)
- Test error conditions and edge cases
- Validate data transformations with known inputs/outputs

### Integration Testing
- Test end-to-end processing with sample VCF files
- Verify output compatibility with original implementation
- Test with various VCF formats and edge cases
- Performance benchmarking against original implementation

### Property-Based Testing
- Generate random VCF-like data for transformation testing
- Test invariants in genotype conversion
- Validate round-trip properties where applicable

## Correctness Properties

*A property is a characteristic or behavior that should hold true across all valid executions of a system-essentially, a formal statement about what the system should do. Properties serve as the bridge between human-readable specifications and machine-verifiable correctness guarantees.*

### Property 1: Variant ID Generation Consistency
*For any* VCF variant record, the generated variant ID should follow the format "CHROM_POS" and be deterministic across multiple runs
**Validates: Requirements 1.2**

### Property 2: Output Format Equivalence
*For any* valid VCF file and target ID file, the new implementation should produce output tables with identical content and structure to the original bcftools-based implementation
**Validates: Requirements 1.5, 7.1, 7.2**

### Property 3: Error Handling with Continuation
*For any* VCF file containing invalid records, the processor should log appropriate warnings and continue processing valid records without terminating
**Validates: Requirements 2.3, 6.1**

### Property 4: Logging Behavior
*For any* processing operation, appropriate log messages should be generated at the correct log levels based on the configured verbosity
**Validates: Requirements 2.7**

### Property 5: Batch Size Respect
*For any* configured batch size, the processor should process variants in groups not exceeding the specified batch size
**Validates: Requirements 4.3**

### Property 6: Progress Indication
*For any* long-running operation, progress indicators should be updated at regular intervals to reflect processing status
**Validates: Requirements 4.4**

### Property 7: Configuration File Parsing
*For any* valid configuration file, all specified parameters should be correctly parsed and applied to the processing configuration
**Validates: Requirements 5.2**

### Property 8: Parameter Validation
*For any* set of input parameters, invalid combinations should be rejected with descriptive error messages before processing begins
**Validates: Requirements 5.3**

### Property 9: Verbosity Mode Control
*For any* verbosity setting (verbose/quiet), the amount and detail of output should correspond appropriately to the selected mode
**Validates: Requirements 5.4**

### Property 10: Dry Run Safety
*For any* dry-run operation, no output files should be created or modified, and a preview of operations should be provided
**Validates: Requirements 5.5**

### Property 11: Invalid Target ID Handling
*For any* target ID file containing invalid entries, those entries should be skipped and reported in the processing summary
**Validates: Requirements 6.2**

### Property 12: Transformation Error Context
*For any* variant transformation that fails, the error message should include specific details about the variant and the reason for failure
**Validates: Requirements 6.3**

### Property 13: VCF Format Validation
*For any* input file, invalid VCF format should be detected and reported before processing begins
**Validates: Requirements 6.4**

### Property 14: Processing Summary Generation
*For any* completed processing run, summary statistics should include counts of total, processed, and skipped variants
**Validates: Requirements 6.5**

### Property 15: Output File Naming Consistency
*For any* output file path, the generated files should follow the expected naming pattern (.gt.txt.gz and .seq.txt.gz)
**Validates: Requirements 7.3**

### Property 16: Output Format Support
*For any* compression setting, the processor should generate valid output files in both compressed and uncompressed formats as requested
**Validates: Requirements 7.4**

### Property 17: Variant Transformation Preservation
*For any* variant (REF, ALT combination), the transformation to annotation format should produce identical results to the original implementation
**Validates: Requirements 7.5**

## Error Handling

The system will implement comprehensive error handling:

1. **VCF File Errors**: Invalid format, corrupted files, missing samples
2. **Target ID Errors**: Invalid format, non-existent IDs
3. **Transformation Errors**: Complex variants that cannot be transformed
4. **I/O Errors**: File permission issues, disk space problems
5. **Memory Errors**: Large file handling with graceful degradation

Each error type will have specific handling strategies:
- Log detailed error information
- Continue processing when possible
- Provide recovery suggestions
- Generate error reports

## Testing Strategy

### Unit Testing
- Test each component in isolation with specific examples
- Mock external dependencies (file I/O, cyvcf2 objects)
- Test error conditions and edge cases
- Validate data transformations with known inputs/outputs
- Focus on critical functions like variant transformation and genotype conversion

### Property-Based Testing
- Use hypothesis or similar library to generate random VCF-like data
- Test variant transformation properties across many input combinations
- Validate genotype conversion invariants
- Test error handling with randomly generated invalid inputs
- Minimum 100 iterations per property test
- Each property test tagged with: **Feature: vcf-processor-optimization, Property {number}: {property_text}**

### Integration Testing
- Test end-to-end processing with sample VCF files
- Verify output compatibility with original implementation
- Test with various VCF formats and edge cases
- Performance benchmarking against original implementation

### Dual Testing Approach
Both unit tests and property tests are necessary and complementary:
- Unit tests verify specific examples and edge cases work correctly
- Property tests verify universal properties hold across all inputs
- Together they provide comprehensive coverage of both concrete behavior and general correctness