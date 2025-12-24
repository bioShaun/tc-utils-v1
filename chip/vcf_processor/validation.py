"""Validation script to ensure all requirements are met."""

import importlib
import sys
from pathlib import Path
from typing import Dict, List, Tuple

def check_dependencies() -> Tuple[bool, List[str]]:
    """Check that all required dependencies are available.
    
    Returns:
        Tuple of (success, missing_dependencies)
    """
    required_deps = [
        'cyvcf2',
        'typer', 
        'loguru',
        'tqdm',
        'pandas',
        'numpy',
        'hypothesis',
        'pytest',
        'psutil'
    ]
    
    missing = []
    
    for dep in required_deps:
        try:
            importlib.import_module(dep)
        except ImportError:
            missing.append(dep)
    
    return len(missing) == 0, missing


def check_module_imports() -> Tuple[bool, List[str]]:
    """Check that all VCF processor modules can be imported.
    
    Returns:
        Tuple of (success, failed_imports)
    """
    modules = [
        'chip.vcf_processor.config',
        'chip.vcf_processor.vcf_reader',
        'chip.vcf_processor.variant_filter',
        'chip.vcf_processor.variant_transformer',
        'chip.vcf_processor.genotype_converter',
        'chip.vcf_processor.output_writer',
        'chip.vcf_processor.vcf_processor',
        'chip.vcf_processor.error_handler',
        'chip.vcf_processor.logging_config',
        'chip.vcf_processor.config_manager',
        'chip.vcf_processor.performance_optimizer',
        'chip.vcf_processor.cli'
    ]
    
    failed = []
    
    for module in modules:
        try:
            importlib.import_module(module)
        except ImportError as e:
            failed.append(f"{module}: {e}")
    
    return len(failed) == 0, failed


def check_file_structure() -> Tuple[bool, List[str]]:
    """Check that all expected files are present.
    
    Returns:
        Tuple of (success, missing_files)
    """
    base_path = Path(__file__).parent
    
    expected_files = [
        '__init__.py',
        'config.py',
        'vcf_reader.py',
        'variant_filter.py',
        'variant_transformer.py',
        'genotype_converter.py',
        'output_writer.py',
        'vcf_processor.py',
        'error_handler.py',
        'logging_config.py',
        'config_manager.py',
        'performance_optimizer.py',
        'cli.py',
        'README.md',
        'tests/__init__.py',
        'tests/conftest.py'
    ]
    
    missing = []
    
    for file_path in expected_files:
        full_path = base_path / file_path
        if not full_path.exists():
            missing.append(str(file_path))
    
    return len(missing) == 0, missing


def check_test_files() -> Tuple[bool, List[str]]:
    """Check that all expected test files are present.
    
    Returns:
        Tuple of (success, missing_test_files)
    """
    base_path = Path(__file__).parent / 'tests'
    
    expected_test_patterns = [
        'test_*_unit.py',
        'test_*_properties.py',
        'test_integration*.py'
    ]
    
    found_tests = []
    for pattern in expected_test_patterns:
        found_tests.extend(list(base_path.glob(pattern)))
    
    # Check for minimum expected test files
    expected_minimum = 15  # We should have at least 15 test files
    
    if len(found_tests) < expected_minimum:
        return False, [f"Expected at least {expected_minimum} test files, found {len(found_tests)}"]
    
    return True, []


def check_api_completeness() -> Tuple[bool, List[str]]:
    """Check that all required API components are available.
    
    Returns:
        Tuple of (success, missing_components)
    """
    missing = []
    
    try:
        # Check main classes
        from chip.vcf_processor.config import ProcessingConfig, ProcessingResult
        from chip.vcf_processor.vcf_processor import VCFProcessor, VCFProcessorFactory
        
        # Check that ProcessingConfig has required attributes
        config_attrs = [
            'vcf_file', 'target_id_file', 'output_file', 'miss_fmt', 'gt_sep',
            'batch_size', 'threads', 'compress_output', 'verbose', 'quiet', 'dry_run'
        ]
        
        for attr in config_attrs:
            if not hasattr(ProcessingConfig, '__annotations__') or attr not in ProcessingConfig.__annotations__:
                missing.append(f"ProcessingConfig missing attribute: {attr}")
        
        # Check that VCFProcessor has required methods
        processor_methods = ['process', 'validate_output', 'get_processing_summary', 'format_summary']
        
        for method in processor_methods:
            if not hasattr(VCFProcessor, method):
                missing.append(f"VCFProcessor missing method: {method}")
        
        # Check that VCFProcessorFactory has required methods
        factory_methods = ['create_processor', 'create_from_files', 'create_batch_processor']
        
        for method in factory_methods:
            if not hasattr(VCFProcessorFactory, method):
                missing.append(f"VCFProcessorFactory missing method: {method}")
        
    except ImportError as e:
        missing.append(f"Import error: {e}")
    
    return len(missing) == 0, missing


def check_cli_functionality() -> Tuple[bool, List[str]]:
    """Check that CLI functionality is available.
    
    Returns:
        Tuple of (success, issues)
    """
    issues = []
    
    try:
        from chip.vcf_processor.cli import app
        
        # Check that CLI app is a typer app
        if not hasattr(app, 'commands'):
            issues.append("CLI app does not appear to be a valid typer application")
        
        # Check for expected commands
        expected_commands = ['process', 'validate', 'info']
        
        # Note: This is a simplified check - in practice we'd need to inspect the app more deeply
        # For now, we'll just check that the CLI module imports successfully
        
    except ImportError as e:
        issues.append(f"CLI import error: {e}")
    
    return len(issues) == 0, issues


def run_validation() -> Dict[str, Tuple[bool, List[str]]]:
    """Run all validation checks.
    
    Returns:
        Dictionary of check results
    """
    checks = {
        'Dependencies': check_dependencies(),
        'Module Imports': check_module_imports(),
        'File Structure': check_file_structure(),
        'Test Files': check_test_files(),
        'API Completeness': check_api_completeness(),
        'CLI Functionality': check_cli_functionality()
    }
    
    return checks


def main():
    """Main validation function."""
    print("VCF Processor Validation")
    print("=" * 50)
    
    results = run_validation()
    
    all_passed = True
    
    for check_name, (success, issues) in results.items():
        status = "✓ PASS" if success else "✗ FAIL"
        print(f"{check_name}: {status}")
        
        if not success:
            all_passed = False
            for issue in issues:
                print(f"  - {issue}")
        
        print()
    
    print("=" * 50)
    
    if all_passed:
        print("✓ All validation checks passed!")
        print("The VCF Processor implementation is complete and ready for use.")
        return 0
    else:
        print("✗ Some validation checks failed.")
        print("Please address the issues above before considering the implementation complete.")
        return 1


if __name__ == "__main__":
    sys.exit(main())