# Plan: Refactor VCF Genotype Statistics Script

## Phase 1: Setup and Relocation
- [x] Task: Relocate script from `draft/vcf_genotype_stats.py` to `vcf/vcf_genotype_stats.py` and ensure `vcf/__init__.py` exists. 5ef42e7
- [x] Task: Create initial test directory and file `vcf/tests/test_vcf_genotype_stats.py`. 896d410
- [ ] Task: Conductor - User Manual Verification 'Phase 1: Setup and Relocation' (Protocol in workflow.md)

## Phase 2: Core Logic Refactoring (TDD)
- [ ] Task: **Write Tests** for genotype classification logic (hom_ref, het, hom_alt, missing).
- [ ] Task: **Implement** refactored genotype classification logic using `loguru`.
- [ ] Task: **Write Tests** for VCF processing and aggregation (using a small sample VCF).
- [ ] Task: **Implement** refactored VCF processing logic with `typer` CLI.
- [ ] Task: **Write Tests** for CSV output functionality.
- [ ] Task: **Implement** CSV output logic.
- [ ] Task: Conductor - User Manual Verification 'Phase 2: Core Logic Refactoring (TDD)' (Protocol in workflow.md)

## Phase 3: Quality Assurance and Finalization
- [ ] Task: Verify >80% test coverage for `vcf/vcf_genotype_stats.py`.
- [ ] Task: Run `ruff` and `pyright` checks and resolve any findings.
- [ ] Task: Final manual verification of CLI output formatting and CSV content.
- [ ] Task: Conductor - User Manual Verification 'Phase 3: Quality Assurance and Finalization' (Protocol in workflow.md)
