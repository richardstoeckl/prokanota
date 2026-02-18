#!/usr/bin/env python3
"""
Shared logging utilities for all search and parsing scripts.
All output goes to stdout so Snakemake's shell redirection (2>&1 | tee) captures it.
"""

import logging
import sys
from datetime import datetime


def setup_logger(script_name: str) -> logging.Logger:
    """
    Configure and return a logger that outputs to stdout.
    
    Args:
        script_name: Name of the script using the logger (for log messages)
        
    Returns:
        Configured logger instance
    """
    logger = logging.getLogger(script_name)
    logger.setLevel(logging.INFO)
    
    # Remove any existing handlers to avoid duplicates
    logger.handlers.clear()
    
    # Create handler that outputs to stdout
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    
    # Create formatter with timestamp, level, and message
    formatter = logging.Formatter(
        fmt='%(asctime)s | %(levelname)-8s | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    handler.setFormatter(formatter)
    
    logger.addHandler(handler)
    
    # Prevent propagation to root logger
    logger.propagate = False
    
    return logger


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
