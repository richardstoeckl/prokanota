#!/usr/bin/env python3
"""
MMseqs2 search wrapper (BLAST-like outfmt 6).

Usage example:
    python mmseqs2_search.py \
        --db /path/to/mmseqs2_db \
        --faa proteins.faa \
        --output results.tsv \
        --toolversion tool_version.txt \
        --threads 4 \
        --evalue 1e-3 \
        --sensitivity 5.7 \
        --max-seqs 300
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from datetime import datetime
from logging_utils import setup_logger, log_file_paths, log_parameters, log_statistics, log_command, build_part, log_prokanota_version


def count_fasta_sequences(fasta_path: str) -> int:
    """Count the number of sequences in a FASTA file."""
    count = 0
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


def count_tabular_hits(tsv_path: str) -> int:
    """Count the number of hits in a tabular output file."""
    if not os.path.isfile(tsv_path) or os.path.getsize(tsv_path) == 0:
        return 0
    count = 0
    with open(tsv_path, 'r') as f:
        count = sum(1 for _ in f)
    return count


def main():
    parser = argparse.ArgumentParser(description="MMseqs2 search wrapper (BLAST-like outfmt 6)")
    parser.add_argument("--db", required=True, help="MMseqs2 database path")
    parser.add_argument("--db-name", help="Database name used for logging context")
    parser.add_argument("--faa", required=True, help="FASTA file with protein sequences")
    parser.add_argument("--output", required=True, help="Output file for tabular results")
    parser.add_argument("--toolversion", required=True, help="Output file to record tool version")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use")
    parser.add_argument("--evalue", type=float, default=0.001, help="E-value cutoff")
    parser.add_argument("--sensitivity", type=float, default=None, help="MMseqs2 sensitivity (-s)")
    parser.add_argument("--max-seqs", type=int, default=None, help="Maximum number of target sequences")
    parser.add_argument("--min-seq-id", type=float, default=None, help="Minimum sequence identity (0-1)")
    parser.add_argument("--cov-mode", type=int, default=None, help="Coverage mode (MMseqs2)")
    parser.add_argument("--coverage", type=float, default=None, help="Coverage threshold (0-1)")
    args = parser.parse_args()

    db_name = args.db_name or Path(args.db).stem
    logger = setup_logger(build_part("SEARCH", db_name))
    log_prokanota_version(logger)
    start_time = datetime.now()
    logger.info("=" * 60)
    logger.info("MMseqs2 Search Started")
    logger.info("=" * 60)

    # Log input files
    log_file_paths(logger, database=args.db, input_fasta=args.faa, output_file=args.output)
    
    # Log parameters
    params = {
        'threads': args.threads,
        'evalue': args.evalue,
        'sensitivity': args.sensitivity,
        'max_seqs': args.max_seqs,
        'min_seq_id': args.min_seq_id,
        'cov_mode': args.cov_mode,
        'coverage': args.coverage,
    }
    log_parameters(logger, **{k: v for k, v in params.items() if v is not None})

    if not os.path.exists(args.db):
        logger.error(f"MMseqs2 database path {args.db} not found")
        sys.exit(1)
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
            ["mmseqs", "version"],
            capture_output=True,
            text=True,
            check=True,
        )
        with open(args.toolversion, "w") as tv:
            tv.write(version_result.stdout)
        logger.info(f"Tool version: {version_result.stdout.strip()}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error getting mmseqs2 version: {e}")
        logger.error(f"stderr: {e.stderr}")
        sys.exit(1)

    output_path = Path(args.output)
    tmp_dir = output_path.with_suffix(".mmseqs_tmp")
    tmp_dir.mkdir(parents=True, exist_ok=True)

    format_output = "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

    cmd = [
        "mmseqs",
        "easy-search",
        args.faa,
        args.db,
        args.output,
        str(tmp_dir),
        "--format-output",
        format_output,
        "--threads",
        str(args.threads),
        "--e-value",
        str(args.evalue),
    ]

    if args.sensitivity is not None:
        cmd.extend(["-s", str(args.sensitivity)])
    if args.max_seqs is not None:
        cmd.extend(["--max-seqs", str(args.max_seqs)])
    if args.min_seq_id is not None:
        cmd.extend(["--min-seq-id", str(args.min_seq_id)])
    if args.cov_mode is not None:
        cmd.extend(["--cov-mode", str(args.cov_mode)])
    if args.coverage is not None:
        cmd.extend(["-c", str(args.coverage)])

    log_command(logger, cmd)
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running mmseqs2: {e}")
        logger.error(f"stderr: {e.stderr}")
        sys.exit(1)

    # Count output hits
    output_count = count_tabular_hits(args.output)
    log_statistics(logger, output_hits=output_count)
    
    # Check for no hits
    if output_count == 0 and input_count > 0:
        logger.warning(f"No hits found in database for {input_count} input sequences. Database may not match data or parameters need adjustment.")
    
    # Log execution time
    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    logger.info(f"Execution time: {duration:.2f} seconds")
    logger.info("=" * 60)
    logger.info("MMseqs2 Search Completed Successfully")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
