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


def main():
    parser = argparse.ArgumentParser(description="MMseqs2 search wrapper (BLAST-like outfmt 6)")
    parser.add_argument("--db", required=True, help="MMseqs2 database path")
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

    if not os.path.exists(args.db):
        sys.exit(f"Error: MMseqs2 database path {args.db} not found")
    if not os.path.isfile(args.faa):
        sys.exit(f"Error: FASTA file {args.faa} does not exist")
    if os.path.getsize(args.faa) == 0:
        open(args.output, "w").close()
        sys.exit(0)

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
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Error getting mmseqs2 version: {e}\n")
        sys.stderr.write(f"stderr: {e.stderr}\n")
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

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Error running mmseqs2: {e}\n")
        sys.stderr.write(f"stderr: {e.stderr}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
