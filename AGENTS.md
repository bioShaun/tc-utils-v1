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
