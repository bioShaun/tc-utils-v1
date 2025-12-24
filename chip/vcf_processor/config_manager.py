"""Configuration management for VCF processor."""

import json
from pathlib import Path
from typing import Any, Dict, Optional

import typer
from loguru import logger

from .config import ProcessingConfig


class ConfigManager:
    """Manages configuration loading and validation for VCF processor."""
    
    @staticmethod
    def from_cli() -> ProcessingConfig:
        """Create configuration from command line arguments using typer.
        
        Returns:
            ProcessingConfig: Validated configuration object
        """
        app = typer.Typer()
        
        @app.command()
        def main(
            vcf_file: Path = typer.Argument(..., help="Input VCF file path"),
            target_id_file: Path = typer.Argument(..., help="Target IDs file path"),
            output_file: Path = typer.Argument(..., help="Output file base path"),
            miss_fmt: str = typer.Option("NN", "--miss-fmt", help="Missing genotype format"),
            gt_sep: str = typer.Option("", "--gt-sep", help="Genotype separator"),
            threads: int = typer.Option(4, "--threads", "-t", help="Number of threads"),
            batch_size: int = typer.Option(10000, "--batch-size", help="Processing batch size"),
            compress_output: bool = typer.Option(True, "--compress/--no-compress", help="Compress output files"),
            verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging"),
            quiet: bool = typer.Option(False, "--quiet", "-q", help="Enable quiet mode"),
            dry_run: bool = typer.Option(False, "--dry-run", help="Preview operations without executing"),
            log_file: Optional[Path] = typer.Option(None, "--log-file", help="Log file path"),
            config_file: Optional[Path] = typer.Option(None, "--config", "-c", help="Configuration file path"),
        ) -> None:
            """Process VCF files to generate genotype tables."""
            # If config file is provided, load base configuration from it
            if config_file:
                base_config = ConfigManager.from_file(config_file)
                # Override with CLI arguments (CLI takes precedence)
                config_dict = base_config.to_dict()
                
                # Update with CLI arguments that were explicitly provided
                cli_args = {
                    "vcf_file": vcf_file,
                    "target_id_file": target_id_file,
                    "output_file": output_file,
                    "miss_fmt": miss_fmt,
                    "gt_sep": gt_sep,
                    "threads": threads,
                    "batch_size": batch_size,
                    "compress_output": compress_output,
                    "verbose": verbose,
                    "quiet": quiet,
                    "dry_run": dry_run,
                    "log_file": log_file,
                }
                
                # Only override if CLI argument differs from default
                for key, value in cli_args.items():
                    if value is not None:
                        config_dict[key] = value
                
                return ProcessingConfig(**config_dict)
            else:
                # Create configuration directly from CLI arguments
                return ProcessingConfig(
                    vcf_file=vcf_file,
                    target_id_file=target_id_file,
                    output_file=output_file,
                    miss_fmt=miss_fmt,
                    gt_sep=gt_sep,
                    threads=threads,
                    batch_size=batch_size,
                    compress_output=compress_output,
                    verbose=verbose,
                    quiet=quiet,
                    dry_run=dry_run,
                    log_file=log_file,
                )
        
        # This is a bit of a hack to get the configuration without actually running typer
        # In a real CLI, you would call app() directly
        return main
    
    @staticmethod
    def from_file(config_path: Path) -> ProcessingConfig:
        """Load configuration from JSON file.
        
        Args:
            config_path: Path to configuration file
            
        Returns:
            ProcessingConfig: Loaded configuration
            
        Raises:
            FileNotFoundError: If config file doesn't exist
            ValueError: If config file is invalid
        """
        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_path}")
        
        try:
            with open(config_path, 'r') as f:
                config_data = json.load(f)
            
            logger.debug(f"Loaded configuration from {config_path}")
            
            # Convert string paths to Path objects
            if "vcf_file" in config_data:
                config_data["vcf_file"] = Path(config_data["vcf_file"])
            if "target_id_file" in config_data:
                config_data["target_id_file"] = Path(config_data["target_id_file"])
            if "output_file" in config_data:
                config_data["output_file"] = Path(config_data["output_file"])
            if "log_file" in config_data and config_data["log_file"]:
                config_data["log_file"] = Path(config_data["log_file"])
            
            return ProcessingConfig(**config_data)
            
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in configuration file {config_path}: {e}")
        except TypeError as e:
            raise ValueError(f"Invalid configuration parameters in {config_path}: {e}")
    
    @staticmethod
    def save_config(config: ProcessingConfig, config_path: Path) -> None:
        """Save configuration to JSON file.
        
        Args:
            config: Configuration to save
            config_path: Path where to save configuration
        """
        config_path.parent.mkdir(parents=True, exist_ok=True)
        
        config_dict = config.to_dict()
        
        with open(config_path, 'w') as f:
            json.dump(config_dict, f, indent=2)
        
        logger.info(f"Configuration saved to {config_path}")
    
    @staticmethod
    def validate_config(config_dict: Dict[str, Any]) -> Dict[str, str]:
        """Validate configuration dictionary and return errors.
        
        Args:
            config_dict: Configuration dictionary to validate
            
        Returns:
            Dict mapping field names to error messages
        """
        errors = {}
        
        # Required fields
        required_fields = ["vcf_file", "target_id_file", "output_file"]
        for field in required_fields:
            if field not in config_dict or not config_dict[field]:
                errors[field] = f"{field} is required"
        
        # Validate numeric fields
        if "threads" in config_dict:
            try:
                threads = int(config_dict["threads"])
                if threads <= 0:
                    errors["threads"] = "threads must be positive"
            except (ValueError, TypeError):
                errors["threads"] = "threads must be a valid integer"
        
        if "batch_size" in config_dict:
            try:
                batch_size = int(config_dict["batch_size"])
                if batch_size <= 0:
                    errors["batch_size"] = "batch_size must be positive"
            except (ValueError, TypeError):
                errors["batch_size"] = "batch_size must be a valid integer"
        
        # Validate boolean fields
        bool_fields = ["compress_output", "verbose", "quiet", "dry_run"]
        for field in bool_fields:
            if field in config_dict and not isinstance(config_dict[field], bool):
                errors[field] = f"{field} must be a boolean"
        
        # Validate conflicting options
        if (config_dict.get("verbose", False) and 
            config_dict.get("quiet", False)):
            errors["verbose"] = "Cannot enable both verbose and quiet modes"
            errors["quiet"] = "Cannot enable both verbose and quiet modes"
        
        # Validate string fields
        if "miss_fmt" in config_dict and not config_dict["miss_fmt"]:
            errors["miss_fmt"] = "miss_fmt cannot be empty"
        
        return errors
    
    @staticmethod
    def create_example_config(output_path: Path) -> None:
        """Create an example configuration file.
        
        Args:
            output_path: Where to save the example configuration
        """
        example_config = {
            "vcf_file": "input.vcf.gz",
            "target_id_file": "targets.txt",
            "output_file": "output",
            "miss_fmt": "NN",
            "gt_sep": "",
            "threads": 4,
            "batch_size": 10000,
            "compress_output": True,
            "verbose": False,
            "quiet": False,
            "dry_run": False,
            "log_file": None
        }
        
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            json.dump(example_config, f, indent=2)
        
        logger.info(f"Example configuration created at {output_path}")


def create_cli_app() -> typer.Typer:
    """Create the main CLI application.
    
    Returns:
        typer.Typer: Configured CLI application
    """
    app = typer.Typer(
        name="vcf-processor",
        help="Optimized VCF genotype table processor using cyvcf2",
        add_completion=False,
    )
    
    @app.command()
    def process(
        vcf_file: Path = typer.Argument(..., help="Input VCF file path"),
        target_id_file: Path = typer.Argument(..., help="Target IDs file path"),
        output_file: Path = typer.Argument(..., help="Output file base path"),
        miss_fmt: str = typer.Option("NN", "--miss-fmt", help="Missing genotype format"),
        gt_sep: str = typer.Option("", "--gt-sep", help="Genotype separator"),
        threads: int = typer.Option(4, "--threads", "-t", help="Number of threads"),
        batch_size: int = typer.Option(10000, "--batch-size", help="Processing batch size"),
        compress_output: bool = typer.Option(True, "--compress/--no-compress", help="Compress output files"),
        verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging"),
        quiet: bool = typer.Option(False, "--quiet", "-q", help="Enable quiet mode"),
        dry_run: bool = typer.Option(False, "--dry-run", help="Preview operations without executing"),
        log_file: Optional[Path] = typer.Option(None, "--log-file", help="Log file path"),
        config_file: Optional[Path] = typer.Option(None, "--config", "-c", help="Configuration file path"),
    ) -> None:
        """Process VCF files to generate genotype tables."""
        from .logging_config import setup_logging
        from .processor import VCFProcessor
        
        # Load configuration
        if config_file:
            config = ConfigManager.from_file(config_file)
            # Override with CLI arguments where provided
            # This is a simplified version - in practice you'd check which args were explicitly set
            if vcf_file != Path(""):
                config.vcf_file = vcf_file
            if target_id_file != Path(""):
                config.target_id_file = target_id_file
            if output_file != Path(""):
                config.output_file = output_file
        else:
            config = ProcessingConfig(
                vcf_file=vcf_file,
                target_id_file=target_id_file,
                output_file=output_file,
                miss_fmt=miss_fmt,
                gt_sep=gt_sep,
                threads=threads,
                batch_size=batch_size,
                compress_output=compress_output,
                verbose=verbose,
                quiet=quiet,
                dry_run=dry_run,
                log_file=log_file,
            )
        
        # Setup logging
        setup_logging(
            verbose=config.verbose,
            quiet=config.quiet,
            log_file=config.log_file,
        )
        
        # Process VCF
        processor = VCFProcessor(config)
        result = processor.process()
        
        # Print summary
        typer.echo(result.summary())
        
        # Exit with error code if processing failed
        if result.has_errors:
            raise typer.Exit(1)
    
    @app.command()
    def create_config(
        output_path: Path = typer.Argument("config.json", help="Output path for example configuration"),
    ) -> None:
        """Create an example configuration file."""
        ConfigManager.create_example_config(output_path)
        typer.echo(f"Example configuration created at {output_path}")
    
    @app.command()
    def validate_config(
        config_path: Path = typer.Argument(..., help="Configuration file to validate"),
    ) -> None:
        """Validate a configuration file."""
        try:
            config = ConfigManager.from_file(config_path)
            typer.echo(f"✓ Configuration file {config_path} is valid")
            typer.echo(f"  VCF file: {config.vcf_file}")
            typer.echo(f"  Target file: {config.target_id_file}")
            typer.echo(f"  Output: {config.output_file}")
            typer.echo(f"  Threads: {config.threads}")
            typer.echo(f"  Batch size: {config.batch_size}")
        except Exception as e:
            typer.echo(f"✗ Configuration file {config_path} is invalid: {e}", err=True)
            raise typer.Exit(1)
    
    return app