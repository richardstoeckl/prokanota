#!/usr/bin/env python3
"""
RPS-BLAST search wrapper.

Usage example:
    python rpsblast_search.py \
        --db /path/to/cdd.db \
        --faa proteins.faa \
        --output results.tsv \
        --toolversion tool_version.txt \
        --threads 4
"""

import argparse
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from logging_utils import (
    build_part,
    log_command,
    log_file_paths,
    log_parameters,
    log_prokanota_version,
    log_statistics,
    setup_logger,
)


def count_fasta_sequences(fasta_path: str) -> int:
    """Count the number of sequences in a FASTA file."""
    count = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def count_tabular_hits(tsv_path: str) -> int:
    """Count the number of hits in a tabular output file."""
    if not os.path.isfile(tsv_path) or os.path.getsize(tsv_path) == 0:
        return 0
    count = 0
    with open(tsv_path) as f:
        count = sum(1 for _ in f)
    return count


def main():
    parser = argparse.ArgumentParser(description="RPS-BLAST search wrapper")
    parser.add_argument("--db", required=True, help="RPS-BLAST database path")
    parser.add_argument("--db-name", help="Database name used for logging context")
    parser.add_argument(
        "--faa", required=True, help="FASTA file with protein sequences"
    )
    parser.add_argument(
        "--output", required=True, help="Output file for tabular results"
    )
    parser.add_argument(
        "--toolversion", required=True, help="Output file to record tool version"
    )
    parser.add_argument(
        "--threads", type=int, default=1, help="Number of threads to use"
    )
    parser.add_argument("--evalue", type=float, default=0.001, help="E-value cutoff")
    args = parser.parse_args()

    db_name = args.db_name or Path(args.db).stem
    logger = setup_logger(build_part("SEARCH", db_name))
    log_prokanota_version(logger)
    start_time = datetime.now()
    logger.info("=" * 60)
    logger.info("RPS-BLAST Search Started")
    logger.info("=" * 60)

    # Log input files
    log_file_paths(
        logger, database=args.db, input_fasta=args.faa, output_file=args.output
    )

    # Log parameters
    log_parameters(logger, threads=args.threads, evalue=args.evalue)

    # TODO reinstate proper check
    # if not os.path.isfile(args.db + ".pal"):
    #     logger.error(f"RPS-BLAST database {args.db} not found")
    #     sys.exit(1)
    if not os.path.isfile(args.faa):
        logger.error(f"FASTA file {args.faa} does not exist")
        sys.exit(1)

    # Handle empty input
    if os.path.getsize(args.faa) == 0:
        logger.warning("Input FASTA file is empty")
        open(args.output, "w").close()
        logger.info("Created empty output file")
        sys.exit(0)

    # Count input sequences
    input_count = count_fasta_sequences(args.faa)
    logger.info(f"Input sequences: {input_count}")

    # Record tool version
    try:
        version_result = subprocess.run(
            ["rpsblast", "-version"], capture_output=True, text=True, check=True
        )
        with open(args.toolversion, "w") as tv:
            tv.write(version_result.stdout)
        version_text = " | ".join(version_result.stdout.strip().split("\n"))
        logger.info(f"Tool version: {version_text}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error getting rpsblast version: {e}")
        sys.exit(1)

    # Run rpsblast
    cmd = [
        "rpsblast",
        "-query",
        args.faa,
        "-db",
        args.db,
        "-out",
        args.output,
        "-seg",
        "no",
        "-comp_based_stats",
        "1",
        "-outfmt",
        "6",
        "-evalue",
        str(args.evalue),
        "-num_threads",
        str(args.threads),
    ]

    log_command(logger, cmd)

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running rpsblast: {e}")
        logger.error(f"stderr: {e.stderr}")
        sys.exit(1)

    # Count output hits
    output_count = count_tabular_hits(args.output)
    log_statistics(logger, output_hits=output_count)

    # Check for no hits
    if output_count == 0 and input_count > 0:
        logger.warning(
            f"No hits found in database for {input_count} input sequences. Database may not match data or parameters need adjustment."
        )

    # Log execution time
    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    logger.info(f"Execution time: {duration:.2f} seconds")
    logger.info("=" * 60)
    logger.info("RPS-BLAST Search Completed Successfully")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
