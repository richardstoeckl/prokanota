#!/usr/bin/env python3
"""
Shared logging utilities for all search and parsing scripts.
All output goes to stdout so Snakemake's shell redirection (2>&1 | tee) captures it.
"""

import logging
import re
import sys
from pathlib import Path


def get_prokanota_version() -> str:
    """
    Read version from prokanota.VERSION file.
    
    Returns the version string, or "unknown" if the file cannot be read.
    This function enables scripts to run in isolated conda environments
    where the prokanota package may not be installed.
    """
    try:
        version_file = Path(__file__).parent.parent.parent / "prokanota.VERSION"
        return version_file.read_text(encoding="utf-8").strip()
    except (FileNotFoundError, IOError):
        return "unknown"


def normalize_part_token(value: str) -> str:
    """Normalize a part token to uppercase snake-style text for stable parsing."""
    token = re.sub(r"[^A-Za-z0-9]+", "_", value.strip()).strip("_")
    return token.upper() if token else "UNKNOWN"


def build_part(stage: str, name: str | None = None) -> str:
    """Build a stable Part field, e.g. SEARCH_ARCOG or FINALIZE_MERGE."""
    stage_token = normalize_part_token(stage)
    if not name:
        return stage_token
    return f"{stage_token}_{normalize_part_token(name)}"


def setup_logger(part_name: str, level: int = logging.INFO) -> logging.Logger:
    """
    Configure and return a logger that outputs to stdout.
    
    Args:
        part_name: Part field value used by downstream parsing
        level: Logging level for both the logger and its stdout handler
        
    Returns:
        Configured logger instance
    """
    logger = logging.getLogger(part_name)
    logger.setLevel(level)
    
    # Remove any existing handlers to avoid duplicates
    logger.handlers.clear()
    
    # Create handler that outputs to stdout
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(level)
    
    # Create formatter with timestamp, part, level, and message
    formatter = logging.Formatter(
        fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    handler.setFormatter(formatter)
    
    logger.addHandler(handler)
    
    # Prevent propagation to root logger
    logger.propagate = False

    return logger


def log_prokanota_version(logger: logging.Logger) -> None:
    """Log the prokanota version with the current logger Part context."""
    logger.info(f"prokanota version: {get_prokanota_version()}")


def log_command(logger: logging.Logger, command: list) -> None:
    """
    Log a command that will be executed.
    
    Args:
        logger: Logger instance
        command: List of command arguments
    """
    logger.info(f"Executing: {' '.join(str(c) for c in command)}")


def log_file_paths(logger: logging.Logger, **paths) -> None:
    """
    Log file paths being used.
    
    Args:
        logger: Logger instance
        **paths: Keyword arguments like input_file="/path/to/file", db="/path/to/db"
    """
    for name, path in paths.items():
        logger.info(f"{name}: {path}")


def log_parameters(logger: logging.Logger, **params) -> None:
    """
    Log search/parsing parameters.
    
    Args:
        logger: Logger instance
        **params: Keyword arguments like evalue=0.001, threads=4
    """
    for name, value in params.items():
        logger.info(f"Parameter: {name}={value}")


def log_statistics(logger: logging.Logger, **stats) -> None:
    """
    Log execution statistics.
    
    Args:
        logger: Logger instance
        **stats: Keyword arguments like input_sequences=1000, output_hits=250
    """
    for name, value in stats.items():
        logger.info(f"Statistic: {name}={value}")
