#!/usr/bin/env python3
"""
DIAMOND search wrapper.

Usage example:
    python diamond_search.py \
        --db /path/to/database.dmnd \
        --faa proteins.faa \
        --output results.tsv \
        --toolversion tool_version.txt \
        --threads 4
"""

import argparse
import subprocess
import sys
import os


def main():
    parser = argparse.ArgumentParser(description="DIAMOND search wrapper")
    parser.add_argument("--db", required=True, help="DIAMOND database path")
    parser.add_argument("--faa", required=True, help="FASTA file with protein sequences")
    parser.add_argument("--output", required=True, help="Output file for tabular results")
    parser.add_argument("--toolversion", required=True, help="Output file to record tool version")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use")
    parser.add_argument("--evalue", type=float, default=0.001, help="E-value cutoff")
    args = parser.parse_args()

    if not os.path.isfile(args.db):
        sys.exit(f"Error: DIAMOND database {args.db} not found")
    if not os.path.isfile(args.faa):
        sys.exit(f"Error: FASTA file {args.faa} does not exist")
    if os.path.getsize(args.faa) == 0:
        open(args.output, "w").close()
        sys.exit(0)

    # Record tool version
    try:
        version_result = subprocess.run(
            ["diamond", "version"],
            capture_output=True,
            text=True,
            check=True
        )
        with open(args.toolversion, "w") as tv:
            tv.write(version_result.stdout)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Error getting diamond version: {e}\n")
        sys.exit(1)

    # Run diamond blastp
    cmd = [
        "diamond",
        "blastp",
        "--query", args.faa,
        "--db", args.db,
        "--out", args.output,
        "--outfmt", "6",
        "--evalue", str(args.evalue),
        "--threads", str(args.threads)
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Error running diamond: {e}\n")
        sys.stderr.write(f"stderr: {e.stderr}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
