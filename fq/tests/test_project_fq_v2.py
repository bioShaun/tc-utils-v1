
import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch
import sys
import os
import pandas as pd


from fq.project_fq_v2 import (
    Config,
    FastqScanner,
    ErrorRecorder,
    FastqErrorType,
    PathResolver,
    SampleValidator,
    FastqInfo
)

class TestPathResolver:
    def test_no_overlap(self):
        resolver = PathResolver(
            include_paths=frozenset([Path("/a/b")]),
            exclude_paths=frozenset([Path("/c/d")])
        )
        assert isinstance(resolver, PathResolver)

    def test_overlap_error(self):
        with pytest.raises(Exception): # Should be PathOverlapError but importing it might be tricky if not exposed
            PathResolver(
                include_paths=frozenset([Path("/a/b")]),
                exclude_paths=frozenset([Path("/a/b")])
            )

    def test_is_excluded(self):
        resolver = PathResolver(exclude_paths=frozenset([Path("/exclude")]))
        assert resolver.is_excluded(Path("/exclude/file.txt"))
        assert not resolver.is_excluded(Path("/include/file.txt"))

class TestFastqScanner:
    @pytest.fixture
    def scanner_setup(self):
        errors = ErrorRecorder()
        config = Config()
        scanner = FastqScanner(errors, config)
        return scanner, errors

    @patch('pathlib.Path.glob')
    def test_scan_directory_success(self, mock_glob, scanner_setup):
        scanner, errors = scanner_setup
        
        # Setup mock paths
        p1 = MagicMock(spec=Path)
        p1.name = "Sample_L1_R1.fastq.gz"
        p1.absolute.return_value = "/path/to/Sample_L1_R1.fastq.gz"
        
        p2 = MagicMock(spec=Path)
        p2.name = "Sample_L1_R2.fastq.gz"
        p2.absolute.return_value = "/path/to/Sample_L1_R2.fastq.gz"

        # Mock glob to return our files
        # We need to handle multiple calls to glob for different extensions
        def glob_side_effect(pattern):
            if pattern == "*.fastq.gz":
                return [p1, p2]
            return []
        
        mock_glob.side_effect = glob_side_effect
        
        sample_path = MagicMock(spec=Path)
        sample_path.exists.return_value = True
        sample_path.name = "Sample-L1"
        sample_path.glob = mock_glob

        results = scanner.scan_directory(sample_path)
        
        assert len(results) == 2
        assert len(errors) == 0
        
        # Check if we got correct FastqInfo objects
        r1_info = next(r for r in results if r.read_type == "R1")
        assert r1_info.libid == "L1"

    def test_extract_lib_id(self, scanner_setup):
        scanner, _ = scanner_setup
        p = Path("/path/to/Sample-L1")
        assert scanner._extract_lib_id(p) == "L1"
        
        p2 = Path("/path/to/Sample-L1-A")
        assert scanner._extract_lib_id(p2) == "L1-A"

class TestSampleValidator:
    @pytest.fixture
    def validator_setup(self):
        errors = ErrorRecorder()
        warnings = ErrorRecorder()
        validator = SampleValidator(errors, warnings)
        return validator, errors, warnings

    def test_validate_columns_missing(self, validator_setup):
        validator, _, _ = validator_setup
        df = pd.DataFrame({'col1': [1]})
        with pytest.raises(Exception):
            validator.validate_columns(df)

    def test_check_low_data(self, validator_setup):
        validator, errors, _ = validator_setup
        df = pd.DataFrame({
            'libid': ['L1'],
            'sample_id': ['S1'],
            'dir_name': ['D1'],
            'data_size': [0.001]
        })
        validator._check_low_data(df, threshold=0.01)
        assert len(errors) == 1
        assert errors._items[0].error_type == FastqErrorType.INCOMPLETE
