#!/usr/bin/env python3
"""
Run tests for add_flank_sequence.py
"""

import sys
from pathlib import Path

import pytest

# Add the primer directory to the path
primer_dir = Path(__file__).parent.parent
sys.path.insert(0, str(primer_dir))

if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main(
        [str(Path(__file__).parent / "test_add_flank_sequence.py"), "-v", "--tb=short"]
    )
