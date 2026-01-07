# tc-pytools

A collection of Python-based CLI utilities for bioinformatics data processing (specifically GTF/BED file manipulation) and general system automation.

## Project Structure
- `gtf/` - Genomic file manipulation tools
  - `split_bed.py` - Legacy script (Reference only, DO NOT use as template)
  - `split_large_genome_bed_v2.py` - **Gold Standard** (Always reference this!)
- `tests/` - Pytest test suite
- `pyproject.toml` - Dependency management

## Available Skills

Use the `skill` tool to load detailed instructions when needed:

| Skill | Trigger | Usage |
|-------|---------|-------|
| `python-review` | Code review, security analysis, quality assessment | `skill(name="python-review")` |
| `python-refactor` | Refactoring legacy code, applying clean code principles | `skill(name="python-refactor")` |
| `python-modern-cli` | Creating new CLI tools from scratch | `skill(name="python-modern-cli")` |

## External File Loading

CRITICAL: When working on this project, use lazy loading for reference files:
- When refactoring or creating CLI tools, read `@gtf/split_large_genome_bed_v2.py` as the reference implementation
- Do NOT preemptively load all files - only load when the specific task requires it

## Code Standards

### Required Stack
- **CLI Framework:** `typer` with `Annotated` type hints
- **Logging:** `loguru` (NEVER use `print` or standard `logging`)
- **Terminal Output:** `rich.console.Console` for user-facing messages
- **Data Processing:** `polars` (preferred) or `pandas` (legacy/small data only)
- **Paths:** `pathlib.Path` (NEVER use `os.path`)
- **Bioinformatics:** `cyvcf2` (VCF), `pyfaidx` (FASTA), `pysam` (BAM/SAM)

### Bioinformatics Libraries

Use modern, actively maintained packages for bioinformatics data processing:

| Data Type | Primary Package | When to Use | Alternative/Fallback |
|-----------|----------------|-------------|---------------------|
| **VCF files** | `cyvcf2` | All VCF parsing/reading tasks | `pysam.VariantFile` (when cyvcf2 lacks specific functionality) |
| **FASTA files** | `pyfaidx` | Large genomes, random access | `Biopython SeqIO` (small files, format conversion, sequence manipulation) |
| **BAM/SAM files** | `pysam` | All alignment file operations | None |
| **BED/GTF files** | `polars` / `pandas` | Tabular genomic coordinates | Custom parsers (avoid) |

**Version requirements:** See `pyproject.toml` for minimum supported versions.

**VCF Processing:**
```python
# Preferred: cyvcf2 for high-performance reading
from cyvcf2 import VCF

for variant in VCF('variants.vcf.gz'):
    chrom, pos = variant.CHROM, variant.POS
    ref, alt = variant.REF, variant.ALT[0]
    
    # Access genotypes efficiently
    gts = variant.gt_types  # 0=HOM_REF, 1=HET, 2=HOM_ALT, 3=UNKNOWN
    
# Fallback: pysam.VariantFile when cyvcf2 lacks features
from pysam import VariantFile

with VariantFile('output.vcf', 'w', header=custom_header) as vcf_out:
    vcf_out.write(record)
```

**FASTA Processing:**
```python
# Preferred: pyfaidx for large genomes (memory-efficient, indexed access)
from pyfaidx import Fasta

genome = Fasta('genome.fa')
sequence = genome['chr1'][1000:2000]  # O(1) random access
reverse_comp = genome['chr1'][1000:2000].reverse.complement

# Acceptable: Biopython for specific use cases
from Bio import SeqIO
from Bio.Seq import Seq

# ✓ Good: Small files, format conversion
records = list(SeqIO.parse('small.fasta', 'fasta'))
SeqIO.convert('input.fasta', 'fasta', 'output.gb', 'genbank')

# ✓ Good: Sequence manipulation
protein = Seq('ATGGCCATTGTAATG').translate()

# ✗ Avoid: Large genome parsing (use pyfaidx instead)
genome = {rec.id: rec.seq for rec in SeqIO.parse('genome.fa', 'fasta')}  # Memory inefficient!
```

**BAM/SAM Processing:**
```python
import pysam

# Region-based fetching (memory-efficient)
with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for read in bam.fetch('chr1', 1000, 2000):
        if read.mapping_quality >= 30:
            print(f"{read.query_name}: {read.reference_start}")
            
# Index operations
pysam.index('input.bam')
```

**Performance Guidelines:**
- **cyvcf2** is 5-10× faster than PyVCF for large VCF files
- **pyfaidx** uses minimal memory regardless of genome size (via FAIDX indexing)
- **pysam** handles BAM/CRAM efficiently with built-in decompression

**When Biopython is appropriate:**
- File size < 100MB
- Format conversions (FASTA ↔ GenBank, etc.)
- Sequence manipulation (translation, reverse complement, motif finding)
- Phylogenetics and alignment tasks (Phylo, AlignIO modules)

### Type Hints (Mandatory)
```python
# Use Python 3.10+ syntax
def process(items: list[str], config: dict[str, Any] | None = None) -> Path:
    ...
```

### Docstrings (Required)
Use Google or NumPy style for all public functions:
```python
def validate_file(path: Path) -> bool:
    """
    Validate that the file exists and is readable.
    
    Args:
        path: Path to the file to validate
        
    Returns:
        True if file is valid
        
    Raises:
        typer.Exit: If validation fails
    """
```

## Error Handling Patterns

Follow this standard pattern for all CLI tools:

```python
from loguru import logger
from rich.console import Console

console = Console()

def validate_input(path: Path) -> None:
    """Validate input with proper error handling."""
    if not path.exists():
        logger.error(f"File not found: {path}")
        console.print(f"[red]Error:[/red] File not found: {path}")
        raise typer.Exit(code=1)
```

**Rules:**
- Log errors with `logger.error()` for debugging
- Display user-friendly messages with `console.print()` using Rich markup
- Use `[red]Error:[/red]` for errors, `[yellow]Warning:[/yellow]` for warnings, `[green]✓[/green]` for success
- Always `raise typer.Exit(code=1)` on failure, never `sys.exit()`

## Testing Standards

### Framework & Commands
```bash
pytest tests/                    # Run all tests
pytest tests/ -v                 # Verbose output
pytest tests/ --cov=gtf          # With coverage
```

### Naming Conventions
- Test files: `test_<module_name>.py`
- Test functions: `test_<function_name>_<scenario>()`
- Fixtures: Use `conftest.py` for shared fixtures

### Structure
```python
def test_validate_bed_file_returns_error_for_missing_file(tmp_path: Path) -> None:
    """Test that validation fails for non-existent files."""
    # Arrange
    fake_path = tmp_path / "nonexistent.bed"
    
    # Act & Assert
    with pytest.raises(SystemExit):
        validate_bed_file(fake_path)
```

## Git Conventions

### Commit Message Format
```
<type>: <short description>

<optional body>
```

**Types:**
- `feat`: New feature
- `fix`: Bug fix
- `refactor`: Code refactoring (no functional change)
- `docs`: Documentation only
- `test`: Adding or updating tests
- `chore`: Maintenance tasks

**Examples:**
```
feat: add coordinate transformation for split BED files
refactor: modernize split_bed.py with typer and loguru
fix: handle empty BED files gracefully
```

## Development Workflows

### Refactoring Legacy Scripts
When encountering scripts using `argparse`, `os.path`, or `print`:
1. Load the `python-refactor` skill
2. Read `gtf/split_large_genome_bed_v2.py` as reference
3. **Analyze** the existing logic and identify anti-patterns
4. **Modernize** using `typer`, `loguru`, `rich`, and `pathlib`
5. **Validate** all inputs explicitly with proper error messages
6. **Test** the refactored script

### Creating New Tools
1. Load the `python-modern-cli` skill
2. Read `gtf/split_large_genome_bed_v2.py` as reference
3. Start with a `typer` skeleton using `Annotated` arguments
4. Add validation functions early
5. Setup logging in `main()` based on `--verbose` flag
6. Write tests alongside implementation

## Common Issues

| Issue | Solution |
|-------|----------|
| Import errors for `typer`/`pandas` | Check `pyproject.toml` dependencies |
| Processing >1GB files is slow | Switch from `pandas` to `polars` |
| Type hints not working | Ensure Python 3.10+ and proper imports |
| Rich colors not showing | Check terminal supports ANSI colors |
| Biopython SeqIO is slow on large FASTA | Use `pyfaidx` for indexed random access |
| VCF parsing memory errors | Switch to `cyvcf2` with streaming iteration |
| BAM file crashes or hangs | Use `pysam` with region-based `fetch()` instead of iterating all reads |
