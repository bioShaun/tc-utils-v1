"""Logging configuration for VCF processor."""

import logging
import sys
from pathlib import Path
from typing import Optional

from loguru import logger


def setup_logging(
    verbose: bool = False,
    quiet: bool = False,
    log_file: Optional[Path] = None,
) -> None:
    """Set up logging configuration.
    
    Args:
        verbose: Enable verbose logging (DEBUG level)
        quiet: Enable quiet mode (WARNING level only)
        log_file: Optional log file path
    """
    # Remove default logger
    logger.remove()
    
    # Determine log level
    if quiet:
        level = "WARNING"
    elif verbose:
        level = "DEBUG"
    else:
        level = "INFO"
    
    # Console handler
    logger.add(
        sys.stderr,
        level=level,
        format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | "
               "<level>{level: <8}</level> | "
               "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - "
               "<level>{message}</level>",
        colorize=True,
    )
    
    # File handler if specified
    if log_file:
        logger.add(
            log_file,
            level="DEBUG",  # Always log everything to file
            format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | "
                   "{name}:{function}:{line} - {message}",
            rotation="10 MB",
            retention="7 days",
        )
    
    # Suppress some noisy loggers
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("urllib3").setLevel(logging.WARNING)