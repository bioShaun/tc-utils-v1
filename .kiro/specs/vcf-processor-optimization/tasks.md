# Implementation Plan: VCF Processor Optimization

## Overview

This implementation plan converts the existing bcftools-based VCF processor to use cyvcf2, improves code structure and maintainability, and adds comprehensive testing. The approach focuses on maintaining backward compatibility while modernizing the codebase.

## Tasks

- [x] 1. Set up project structure and dependencies
  - Create new module structure with proper package organization
  - Add cyvcf2, pytest, hypothesis, and other required dependencies to requirements
  - Set up logging configuration and basic project scaffolding
  - _Requirements: 2.1, 2.7_

- [x] 2. Implement core data models and configuration
  - [x] 2.1 Create data models using dataclasses
    - Define VariantInfo, ProcessingConfig, and ProcessingResult dataclasses
    - Add type hints and validation methods
    - _Requirements: 2.2, 2.4_

  - [x] 2.2 Write property test for data model validation
    - **Property 8: Parameter Validation**
    - **Validates: Requirements 5.3**

  - [x] 2.3 Implement ConfigManager class
    - Add CLI argument parsing with typer
    - Add configuration file support
    - Add parameter validation
    - _Requirements: 5.1, 5.2, 5.3_

  - [x] 2.4 Write unit tests for ConfigManager
    - Test CLI parsing with various argument combinations
    - Test configuration file loading and validation
    - _Requirements: 5.1, 5.2, 5.3_

- [x] 3. Implement VCF reading with cyvcf2
  - [x] 3.1 Create VCFReader class
    - Implement cyvcf2-based VCF file reading
    - Add sample name extraction
    - Add variant iteration with optional filtering
    - _Requirements: 1.1, 1.4_

  - [x] 3.2 Write property test for VCF reading
    - **Property 1: Variant ID Generation Consistency**
    - **Validates: Requirements 1.2**

  - [x] 3.3 Implement VariantFilter class
    - Add target ID loading from file
    - Add variant filtering logic
    - _Requirements: 1.3, 6.2_

  - [x] 3.4 Write property test for variant filtering
    - **Property 11: Invalid Target ID Handling**
    - **Validates: Requirements 6.2**

- [x] 4. Checkpoint - Ensure VCF reading tests pass
  - Ensure all tests pass, ask the user if questions arise.

- [x] 5. Implement variant transformation logic
  - [x] 5.1 Create VariantTransformer class
    - Port existing transformOneAlt and transformAlt functions
    - Add comprehensive error handling for complex variants
    - Add validation methods
    - _Requirements: 6.3, 7.5_

  - [x] 5.2 Write property test for variant transformation
    - **Property 17: Variant Transformation Preservation**
    - **Validates: Requirements 7.5**

  - [x] 5.3 Write property test for transformation error handling
    - **Property 12: Transformation Error Context**
    - **Validates: Requirements 6.3**

  - [x] 5.4 Write unit tests for variant transformation edge cases
    - Test SNPs, insertions, deletions, and complex variants
    - Test error conditions with invalid inputs
    - _Requirements: 6.3, 7.5_

- [x] 6. Implement genotype conversion
  - [x] 6.1 Create GenotypeConverter class
    - Port existing genotype conversion logic
    - Add batch processing capabilities
    - Optimize with vectorized operations where possible
    - _Requirements: 4.2, 7.1, 7.2_

  - [x] 6.2 Write property test for genotype conversion
    - **Property 2: Output Format Equivalence**
    - **Validates: Requirements 1.5, 7.1, 7.2**

  - [x] 6.3 Implement batch processing with configurable sizes
    - Add streaming processing for large files
    - Add progress tracking
    - _Requirements: 4.1, 4.3, 4.4_

  - [x] 6.4 Write property test for batch processing
    - **Property 5: Batch Size Respect**
    - **Validates: Requirements 4.3**

- [x] 7. Implement output management
  - [x] 7.1 Create OutputWriter class with context manager support
    - Add support for compressed and uncompressed output
    - Implement proper file naming conventions
    - Add append mode for batch processing
    - _Requirements: 2.6, 7.3, 7.4_

  - [x] 7.2 Write property test for output formatting
    - **Property 15: Output File Naming Consistency**
    - **Validates: Requirements 7.3**

  - [x] 7.3 Write property test for compression support
    - **Property 16: Output Format Support**
    - **Validates: Requirements 7.4**

- [-] 8. Checkpoint - Ensure core functionality tests pass
  - Ensure all tests pass, ask the user if questions arise.

- [x] 9. Implement main VCFProcessor orchestrator
  - [x] 9.1 Create VCFProcessor main class
    - Integrate all components into main processing pipeline
    - Add comprehensive error handling and logging
    - Add progress reporting and summary statistics
    - _Requirements: 2.3, 2.7, 6.1, 6.4, 6.5_

  - [x] 9.2 Write property test for error handling
    - **Property 3: Error Handling with Continuation**
    - **Validates: Requirements 2.3, 6.1**

  - [x] 9.3 Write property test for logging behavior
    - **Property 4: Logging Behavior**
    - **Validates: Requirements 2.7**

  - [x] 9.4 Write property test for VCF validation
    - **Property 13: VCF Format Validation**
    - **Validates: Requirements 6.4**

- [ ] 10. Implement CLI interface and modes
  - [ ] 10.1 Add verbose/quiet modes and dry-run support
    - Implement output verbosity controls
    - Add dry-run mode that previews operations without executing
    - _Requirements: 5.4, 5.5_

  - [ ] 10.2 Write property test for verbosity modes
    - **Property 9: Verbosity Mode Control**
    - **Validates: Requirements 5.4**

  - [ ] 10.3 Write property test for dry-run safety
    - **Property 10: Dry Run Safety**
    - **Validates: Requirements 5.5**

  - [ ] 10.4 Add progress indicators and summary reporting
    - Implement progress bars for long operations
    - Generate comprehensive processing summaries
    - _Requirements: 4.4, 6.5_

  - [ ] 10.5 Write property test for progress indication
    - **Property 6: Progress Indication**
    - **Validates: Requirements 4.4**

- [ ] 11. Integration testing and compatibility verification
  - [ ] 11.1 Create integration tests with sample data
    - Test end-to-end processing with various VCF files
    - Compare outputs with original implementation
    - _Requirements: 1.5, 7.1, 7.2, 7.5_

  - [ ] 11.2 Write property test for processing summary
    - **Property 14: Processing Summary Generation**
    - **Validates: Requirements 6.5**

  - [ ] 11.3 Write integration property tests
    - Test complete pipeline with randomly generated valid inputs
    - Verify output consistency across multiple runs
    - _Requirements: 1.5, 7.1, 7.2_

- [ ] 12. Performance optimization and final polish
  - [ ] 12.1 Optimize performance-critical sections
    - Profile code and optimize bottlenecks
    - Ensure memory usage stays reasonable for large files
    - _Requirements: 4.1, 4.5_

  - [ ] 12.2 Add comprehensive documentation and examples
    - Add docstrings following Google style
    - Create usage examples and migration guide
    - _Requirements: 2.5_

  - [ ] 12.3 Final code review and cleanup
    - Ensure Google Python style compliance
    - Remove any remaining bcftools dependencies
    - Verify all requirements are met
    - _Requirements: 1.1, 2.5_

- [ ] 13. Final checkpoint - Ensure all tests pass
  - Ensure all tests pass, ask the user if questions arise.

## Notes

- Each task references specific requirements for traceability
- Checkpoints ensure incremental validation
- Property tests validate universal correctness properties
- Unit tests validate specific examples and edge cases
- The implementation maintains backward compatibility while modernizing the codebase