# Spec: Refactor VCF Genotype Statistics Script

## Overview
This track involves refactoring the experimental script `draft/vcf_genotype_stats.py` into a production-ready tool within the `vcf/` module. The tool calculates genotype frequencies (hom_ref, het, hom_alt, missing) across all samples in a VCF file.

## Objectives
- Move the script from `draft/vcf_genotype_stats.py` to `vcf/vcf_genotype_stats.py`.
- Align the code with the project's quality standards:
    - Use `typer` for the CLI.
    - Use `loguru` for logging.
    - Implement functional decomposition.
    - Add comprehensive unit tests.
    - Ensure >80% code coverage.
- Maintain existing performance optimizations for `cyvcf2`.

## Functional Requirements
- **Input:** A VCF file (standard or compressed).
- **Genotype Classification:**
    - `hom_ref`: 0/0
    - `het`: 0/1, 1/2, etc.
    - `hom_alt`: 1/1, 2/2, etc.
    - `missing`: Contains `.` (e.g., `./.`, `0/.`).
- **Summary Output:** Print a detailed summary to the console including:
    - Total calls (Samples * Sites).
    - Counts and percentages for each genotype class.
    - Overall Heterozygosity and Missingness rates.
- **Detailed Output (Optional):** If an output path is provided, save site-level statistics to a CSV file.

## Non-Functional Requirements
- **Performance:** Process large VCFs efficiently by leveraging `cyvcf2`'s low-level access to genotype arrays.
- **Robustness:** Gracefully handle malformed VCFs or non-diploid calls (log warnings).

## Acceptance Criteria
- Script moved to `vcf/vcf_genotype_stats.py`.
- `pytest --cov=vcf/vcf_genotype_stats.py` shows >80% coverage.
- Code passes all `ruff` and `pyright` checks.
- CLI output matches or improves upon the original script's readability.
