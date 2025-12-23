# Plan: Refactor VCF Genotype Statistics Script

## Phase 1: Setup and Relocation [checkpoint: 29b5970]
- [x] Task: Relocate script from `draft/vcf_genotype_stats.py` to `vcf/vcf_genotype_stats.py` and ensure `vcf/__init__.py` exists. 5ef42e7
- [x] Task: Create initial test directory and file `vcf/tests/test_vcf_genotype_stats.py`. 896d410
- [x] Task: Conductor - User Manual Verification 'Phase 1: Setup and Relocation' (Protocol in workflow.md) 29b5970

## Phase 2: Core Logic Refactoring (TDD) [checkpoint: 99aff50]
- [x] Task: **Write Tests** for genotype classification logic (hom_ref, het, hom_alt, missing). e04fe1f
- [x] Task: **Implement** refactored genotype classification logic using `loguru`. 0ad8ac5
- [x] Task: **Write Tests** for VCF processing and aggregation (using a small sample VCF). c09a261
- [x] Task: **Implement** refactored VCF processing logic with `typer` CLI. bde0a46
- [x] Task: **Write Tests** for CSV output functionality. 5363d3c
- [x] Task: **Implement** CSV output logic. ca4b830
- [x] Task: Conductor - User Manual Verification 'Phase 2: Core Logic Refactoring (TDD)' (Protocol in workflow.md) 99aff50

## Phase 3: Quality Assurance and Finalization
- [x] Task: Verify >80% test coverage for `vcf/vcf_genotype_stats.py`. 8411895
- [x] Task: Run `ruff` and `pyright` checks and resolve any findings. 8c8efc8
- [x] Task: Final manual verification of CLI output formatting and CSV content. 8c8efc8
- [ ] Task: Conductor - User Manual Verification 'Phase 3: Quality Assurance and Finalization' (Protocol in workflow.md)
