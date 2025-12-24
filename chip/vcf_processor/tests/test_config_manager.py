"""Unit tests for ConfigManager class."""

import json
import tempfile
from pathlib import Path

import pytest

from chip.vcf_processor.config import ProcessingConfig
from chip.vcf_processor.config_manager import ConfigManager


class TestConfigManager:
    """Unit tests for ConfigManager class."""
    
    def test_from_file_valid_config(self, temp_dir: Path) -> None:
        """Test loading valid configuration from file."""
        # Create test files
        vcf_file = temp_dir / "test.vcf"
        vcf_file.write_text("##fileformat=VCFv4.2\n")
        
        target_file = temp_dir / "targets.txt"
        target_file.write_text("chr1_100\n")
        
        # Create config file
        config_data = {
            "vcf_file": str(vcf_file),
            "target_id_file": str(target_file),
            "output_file": str(temp_dir / "output"),
            "miss_fmt": "NN",
            "gt_sep": "",
            "threads": 8,
            "batch_size": 5000,
            "compress_output": False,
            "verbose": True,
            "quiet": False,
            "dry_run": True,
            "log_file": str(temp_dir / "log.txt")
        }
        
        config_file = temp_dir / "config.json"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)
        
        # Load configuration
        config = ConfigManager.from_file(config_file)
        
        # Verify all fields are loaded correctly
        assert config.vcf_file == vcf_file
        assert config.target_id_file == target_file
        assert config.output_file == temp_dir / "output"
        assert config.miss_fmt == "NN"
        assert config.gt_sep == ""
        assert config.threads == 8
        assert config.batch_size == 5000
        assert config.compress_output is False
        assert config.verbose is True
        assert config.quiet is False
        assert config.dry_run is True
        assert config.log_file == temp_dir / "log.txt"
    
    def test_from_file_minimal_config(self, temp_dir: Path) -> None:
        """Test loading minimal configuration with defaults."""
        # Create test files
        vcf_file = temp_dir / "test.vcf"
        vcf_file.write_text("##fileformat=VCFv4.2\n")
        
        target_file = temp_dir / "targets.txt"
        target_file.write_text("chr1_100\n")
        
        # Create minimal config file
        config_data = {
            "vcf_file": str(vcf_file),
            "target_id_file": str(target_file),
            "output_file": str(temp_dir / "output"),
        }
        
        config_file = temp_dir / "config.json"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)
        
        # Load configuration
        config = ConfigManager.from_file(config_file)
        
        # Verify required fields and defaults
        assert config.vcf_file == vcf_file
        assert config.target_id_file == target_file
        assert config.output_file == temp_dir / "output"
        assert config.miss_fmt == "NN"  # default
        assert config.threads == 4  # default
        assert config.batch_size == 10000  # default
        assert config.compress_output is True  # default
    
    def test_from_file_nonexistent(self, temp_dir: Path) -> None:
        """Test loading from nonexistent file raises FileNotFoundError."""
        nonexistent_file = temp_dir / "nonexistent.json"
        
        with pytest.raises(FileNotFoundError, match="Configuration file not found"):
            ConfigManager.from_file(nonexistent_file)
    
    def test_from_file_invalid_json(self, temp_dir: Path) -> None:
        """Test loading invalid JSON raises ValueError."""
        config_file = temp_dir / "invalid.json"
        config_file.write_text("{ invalid json }")
        
        with pytest.raises(ValueError, match="Invalid JSON"):
            ConfigManager.from_file(config_file)
    
    def test_from_file_invalid_parameters(self, temp_dir: Path) -> None:
        """Test loading with invalid parameters raises ValueError."""
        # Create valid files first
        vcf_file = temp_dir / "test.vcf"
        vcf_file.write_text("##fileformat=VCFv4.2\n")
        
        target_file = temp_dir / "targets.txt"
        target_file.write_text("chr1_100\n")
        
        config_data = {
            "vcf_file": str(vcf_file),
            "target_id_file": str(target_file),
            "output_file": str(temp_dir / "output"),
            "threads": -1,  # Invalid
        }
        
        config_file = temp_dir / "config.json"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)
        
        with pytest.raises(ValueError, match="Threads must be positive"):
            ConfigManager.from_file(config_file)
    
    def test_save_config(self, temp_dir: Path) -> None:
        """Test saving configuration to file."""
        # Create test files
        vcf_file = temp_dir / "test.vcf"
        vcf_file.write_text("##fileformat=VCFv4.2\n")
        
        target_file = temp_dir / "targets.txt"
        target_file.write_text("chr1_100\n")
        
        # Create configuration
        config = ProcessingConfig(
            vcf_file=vcf_file,
            target_id_file=target_file,
            output_file=temp_dir / "output",
            threads=8,
            verbose=True,
        )
        
        # Save configuration
        config_file = temp_dir / "saved_config.json"
        ConfigManager.save_config(config, config_file)
        
        # Verify file was created and contains correct data
        assert config_file.exists()
        
        with open(config_file, 'r') as f:
            saved_data = json.load(f)
        
        assert saved_data["vcf_file"] == str(vcf_file)
        assert saved_data["target_id_file"] == str(target_file)
        assert saved_data["output_file"] == str(temp_dir / "output")
        assert saved_data["threads"] == 8
        assert saved_data["verbose"] is True
    
    def test_validate_config_valid(self) -> None:
        """Test validation of valid configuration dictionary."""
        config_dict = {
            "vcf_file": "test.vcf",
            "target_id_file": "targets.txt",
            "output_file": "output",
            "threads": 4,
            "batch_size": 1000,
            "compress_output": True,
            "verbose": False,
            "quiet": False,
            "miss_fmt": "NN",
        }
        
        errors = ConfigManager.validate_config(config_dict)
        assert errors == {}
    
    def test_validate_config_missing_required(self) -> None:
        """Test validation with missing required fields."""
        config_dict = {
            "threads": 4,
        }
        
        errors = ConfigManager.validate_config(config_dict)
        
        assert "vcf_file" in errors
        assert "target_id_file" in errors
        assert "output_file" in errors
        assert "vcf_file is required" in errors["vcf_file"]
    
    def test_validate_config_invalid_types(self) -> None:
        """Test validation with invalid field types."""
        config_dict = {
            "vcf_file": "test.vcf",
            "target_id_file": "targets.txt",
            "output_file": "output",
            "threads": "not_a_number",
            "batch_size": -1,
            "compress_output": "not_a_boolean",
            "verbose": "not_a_boolean",
        }
        
        errors = ConfigManager.validate_config(config_dict)
        
        assert "threads" in errors
        assert "batch_size" in errors
        assert "compress_output" in errors
        assert "verbose" in errors
        assert "must be a valid integer" in errors["threads"]
        assert "must be positive" in errors["batch_size"]
        assert "must be a boolean" in errors["compress_output"]
    
    def test_validate_config_conflicting_options(self) -> None:
        """Test validation with conflicting options."""
        config_dict = {
            "vcf_file": "test.vcf",
            "target_id_file": "targets.txt",
            "output_file": "output",
            "verbose": True,
            "quiet": True,  # Conflicting with verbose
        }
        
        errors = ConfigManager.validate_config(config_dict)
        
        assert "verbose" in errors
        assert "quiet" in errors
        assert "Cannot enable both verbose and quiet modes" in errors["verbose"]
    
    def test_validate_config_empty_miss_fmt(self) -> None:
        """Test validation with empty miss_fmt."""
        config_dict = {
            "vcf_file": "test.vcf",
            "target_id_file": "targets.txt",
            "output_file": "output",
            "miss_fmt": "",
        }
        
        errors = ConfigManager.validate_config(config_dict)
        
        assert "miss_fmt" in errors
        assert "cannot be empty" in errors["miss_fmt"]
    
    def test_create_example_config(self, temp_dir: Path) -> None:
        """Test creating example configuration file."""
        example_file = temp_dir / "example.json"
        
        ConfigManager.create_example_config(example_file)
        
        # Verify file was created
        assert example_file.exists()
        
        # Verify content is valid JSON with expected fields
        with open(example_file, 'r') as f:
            example_data = json.load(f)
        
        required_fields = ["vcf_file", "target_id_file", "output_file"]
        for field in required_fields:
            assert field in example_data
        
        # Verify some default values
        assert example_data["miss_fmt"] == "NN"
        assert example_data["threads"] == 4
        assert example_data["batch_size"] == 10000
        assert example_data["compress_output"] is True
    
    def test_create_example_config_creates_directory(self, temp_dir: Path) -> None:
        """Test that creating example config creates parent directories."""
        nested_dir = temp_dir / "nested" / "directory"
        example_file = nested_dir / "example.json"
        
        # Directory doesn't exist yet
        assert not nested_dir.exists()
        
        ConfigManager.create_example_config(example_file)
        
        # Directory should be created
        assert nested_dir.exists()
        assert example_file.exists()


class TestConfigManagerIntegration:
    """Integration tests for ConfigManager with actual file operations."""
    
    def test_save_and_load_roundtrip(self, temp_dir: Path) -> None:
        """Test saving and loading configuration maintains all data."""
        # Create test files
        vcf_file = temp_dir / "test.vcf"
        vcf_file.write_text("##fileformat=VCFv4.2\n")
        
        target_file = temp_dir / "targets.txt"
        target_file.write_text("chr1_100\n")
        
        # Create original configuration
        original_config = ProcessingConfig(
            vcf_file=vcf_file,
            target_id_file=target_file,
            output_file=temp_dir / "output",
            miss_fmt="--",
            gt_sep="/",
            threads=16,
            batch_size=50000,
            compress_output=False,
            verbose=True,
            quiet=False,
            dry_run=True,
            log_file=temp_dir / "test.log",
        )
        
        # Save configuration
        config_file = temp_dir / "roundtrip.json"
        ConfigManager.save_config(original_config, config_file)
        
        # Load configuration
        loaded_config = ConfigManager.from_file(config_file)
        
        # Verify all fields match
        assert loaded_config.vcf_file == original_config.vcf_file
        assert loaded_config.target_id_file == original_config.target_id_file
        assert loaded_config.output_file == original_config.output_file
        assert loaded_config.miss_fmt == original_config.miss_fmt
        assert loaded_config.gt_sep == original_config.gt_sep
        assert loaded_config.threads == original_config.threads
        assert loaded_config.batch_size == original_config.batch_size
        assert loaded_config.compress_output == original_config.compress_output
        assert loaded_config.verbose == original_config.verbose
        assert loaded_config.quiet == original_config.quiet
        assert loaded_config.dry_run == original_config.dry_run
        assert loaded_config.log_file == original_config.log_file
    
    def test_config_with_null_log_file(self, temp_dir: Path) -> None:
        """Test configuration with null log file."""
        # Create test files
        vcf_file = temp_dir / "test.vcf"
        vcf_file.write_text("##fileformat=VCFv4.2\n")
        
        target_file = temp_dir / "targets.txt"
        target_file.write_text("chr1_100\n")
        
        config_data = {
            "vcf_file": str(vcf_file),
            "target_id_file": str(target_file),
            "output_file": str(temp_dir / "output"),
            "log_file": None
        }
        
        config_file = temp_dir / "config.json"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)
        
        # Should load without error
        config = ConfigManager.from_file(config_file)
        assert config.log_file is None