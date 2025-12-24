"""
VCF Processor Optimization Module

This module provides an optimized VCF genotype table processor that replaces
bcftools with cyvcf2 for better performance and maintainability.
"""

from .config import ProcessingConfig, ProcessingResult
from .config_manager import ConfigManager, create_cli_app
from .vcf_reader import VCFReader
from .variant_filter import VariantFilter
from .genotype_converter import GenotypeConverter
from .variant_transformer import VariantTransformer
from .output_writer import OutputWriter
from .vcf_processor import VCFProcessor

__version__ = "1.0.0"
__all__ = [
    "ProcessingConfig",
    "ProcessingResult", 
    "ConfigManager",
    "create_cli_app",
    "VCFReader",
    "VariantFilter",
    "GenotypeConverter",
    "VariantTransformer",
    "OutputWriter",
    "VCFProcessor",
]