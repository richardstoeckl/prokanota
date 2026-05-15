#!/usr/bin/env python3
"""
pyhmmer-based hmmsearch replacement with performance improvements.

This script loads a binary HMM database and pre-fetches all protein sequences from a FASTA
file into memory before running the search. It uses pyhmmer’s multi-threaded hmmsearch
pipeline to query all HMMs against the in-memory sequence block.

Usage example:
    python pyhmmer_hmmsearch.py \
        --db arcogs.h3m \
        --faa proteins.faa \
        --tblout results_tbl.txt \
        [--allresults results_all.txt] \
        [--domtblout results_dom.txt] \
        --toolversion tool_version.txt \
        --threads 4
"""

import argparse
import pyhmmer
import pyhmmer.easel
import pyhmmer.plan7
import sys
import os
import psutil
from pathlib import Path
from datetime import datetime
from logging_utils import setup_logger, log_file_paths, log_parameters, log_statistics, build_part, log_prokanota_version


def count_hmm_models(db_path: str) -> int:
    """Count the number of HMM models in the database."""
    count = 0
    with pyhmmer.plan7.HMMFile(db_path) as hmm_file:
        for _ in hmm_file:
            count += 1
    return count


def main():
    parser = argparse.ArgumentParser(
        description="pyhmmer-based hmmsearch replacement with performance optimizations")
    parser.add_argument("--db", required=True, help="Binary HMM database file (e.g., .h3m file)")
    parser.add_argument("--db-name", help="Database name used for logging context")
    parser.add_argument("--faa", required=True, help="FASTA file with protein sequences")
    parser.add_argument("--tblout", required=True, help="Output file for table results")
    parser.add_argument("--allresults", help="Output file for full results (optional)")
    parser.add_argument("--domtblout", help="Output file for domain table results (optional)")
    parser.add_argument("--toolversion", required=True, help="Output file to record tool version")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads (CPUs) to use")
    args = parser.parse_args()

    db_name = args.db_name or Path(args.db).stem
    logger = setup_logger(build_part("SEARCH", db_name))
    log_prokanota_version(logger)
    start_time = datetime.now()
    logger.info("=" * 60)
    logger.info("pyhmmer Search Started")
    logger.info("=" * 60)

    # Log input files
    files = {'database': args.db, 'input_fasta': args.faa, 'output_tblout': args.tblout}
    if args.allresults:
        files['output_allresults'] = args.allresults
    if args.domtblout:
        files['output_domtblout'] = args.domtblout
    log_file_paths(logger, **files)
    
    # Log parameters
    log_parameters(logger, threads=args.threads)

    # Record tool version
    with open(args.toolversion, "w") as tv:
        tv.write(f"pyhmmer version: {pyhmmer.__version__}\n")
    logger.info(f"Tool version: pyhmmer {pyhmmer.__version__}")

    # Validate and load HMMs
    if not os.path.isfile(args.db):
        logger.error(f"HMM file {args.db} does not exist")
        sys.exit(1)
    
    logger.info("Loading HMM models...")
    with pyhmmer.plan7.HMMFile(args.db) as hmm_file:
        hmms = list(hmm_file)
    
    hmm_count = len(hmms)
    logger.info(f"HMM models loaded: {hmm_count}")
    log_statistics(logger, hmm_models=hmm_count)

    # Pre-fetch sequences
    if not os.path.isfile(args.faa):
        logger.error(f"FASTA file {args.faa} does not exist")
        sys.exit(1)
    
    file_size = os.path.getsize(args.faa) / (1024 * 1024)  # Size in MB
    available_memory = psutil.virtual_memory().available / (1024 * 1024)  # Available memory in MB
    logger.info(f"FASTA file size: {file_size:.2f} MB, Available memory: {available_memory:.2f} MB")
    
    if file_size > available_memory * 0.8:
        logger.warning(f"FASTA file size ({file_size:.2f} MB) may exceed available memory ({available_memory:.2f} MB)")
    
    logger.info("Loading sequences...")
    with pyhmmer.easel.SequenceFile(args.faa, digital=True) as seq_file:
        sequences = seq_file.read_block()
    
    seq_count = len(sequences)
    logger.info(f"Input sequences loaded: {seq_count}")

    # Run hmmsearch
    logger.info("Running hmmsearch...")
    results = pyhmmer.hmmer.hmmsearch(hmms, sequences, cpus=args.threads)

    # Open output files conditionally
    hit_count = 0
    with open(args.tblout, "w") as tblout:
        # Optional outputs
        allresults_file = open(args.allresults, "w") if args.allresults else None
        domtblout_file = open(args.domtblout, "w") if args.domtblout else None

        # Write tblout header without separator
        tblout.write("# target name        query name           accession    E-value  score  bias\n")

        # Optional headers
        if allresults_file:
            allresults_file.write("Full search results:\n")
        if domtblout_file:
            domtblout_file.write("# target name\tp-value\tscore\n")

        # Process results
        for hmm, result in zip(hmms, results):
            hmm_name = hmm.name  # Keep .sr in query name
            hmm_acc = hmm.accession if hmm.accession else "-"
            
            if allresults_file:
                allresults_file.write(f"HMM: {hmm_name} (accession: {hmm_acc})\n")

            for hit in result:
                hit_count += 1
                hit_name = hit.name

                # Full sequence stats
                full_evalue = hit.evalue
                full_score = hit.score
                full_bias = hit.bias

                # Write to tblout with query accession
                tblout.write(f"{hit_name:<20} {hmm_name:<20} {hmm_acc:<12} {full_evalue:<9.1e} {full_score:<6.1f} {full_bias:<5.1f}\n")

                # Optional outputs
                if allresults_file:
                    hit_acc = hit.accession if hit.accession else "-"
                    allresults_file.write(f"  Hit: {hit_name} (accession: {hit_acc})  E-value: {hit.evalue}  Score: {hit.score}  Bias: {hit.bias}\n")
                if domtblout_file:
                    for dom in hit.domains:
                        domtblout_file.write(f"{hit_name}\t{dom.pvalue}\t{dom.score}\n")

            if allresults_file:
                allresults_file.write("\n")

        # Close optional files
        if allresults_file:
            allresults_file.close()
        if domtblout_file:
            domtblout_file.close()

    # Log execution statistics
    log_statistics(logger, output_hits=hit_count)
    
    # Check for no hits
    if hit_count == 0 and seq_count > 0:
        logger.warning(f"No hits found in HMM database for {seq_count} input sequences. Database may not match data or parameters need adjustment.")
    
    # Log execution time
    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    logger.info(f"Execution time: {duration:.2f} seconds")
    logger.info("=" * 60)
    logger.info("pyhmmer Search Completed Successfully")
    logger.info("=" * 60)

if __name__ == "__main__":
    main()