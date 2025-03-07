#!/usr/bin/env python3
"""
Copyright Richard Stöckl 2025.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE or copy at 
https://www.boost.org/LICENSE_1_0.txt)


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

def main():
    parser = argparse.ArgumentParser(
        description="pyhmmer-based hmmsearch replacement with performance optimizations")
    parser.add_argument("--db", required=True, help="Binary HMM database file (e.g., .h3m file)")
    parser.add_argument("--faa", required=True, help="FASTA file with protein sequences")
    parser.add_argument("--tblout", required=True, help="Output file for table results")
    parser.add_argument("--allresults", help="Output file for full results (optional)")
    parser.add_argument("--domtblout", help="Output file for domain table results (optional)")
    parser.add_argument("--toolversion", required=True, help="Output file to record tool version")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads (CPUs) to use")
    args = parser.parse_args()

    # Record tool version
    with open(args.toolversion, "w") as tv:
        tv.write(f"pyhmmer version: {pyhmmer.__version__}\n")

    # Validate and load HMMs
    if not os.path.isfile(args.db):
        sys.exit(f"Error: HMM file {args.db} does not exist")
    with pyhmmer.plan7.HMMFile(args.db) as hmm_file:
        hmms = list(hmm_file)

    # Pre-fetch sequences
    if not os.path.isfile(args.faa):
        sys.exit(f"Error: FASTA file {args.faa} does not exist")
    file_size = os.path.getsize(args.faa) / (1024 * 1024)  # Size in MB
    available_memory = psutil.virtual_memory().available / (1024 * 1024)  # Available memory in MB
    if file_size > available_memory * 0.8:
        sys.stderr.write(f"Warning: FASTA file size ({file_size:.2f} MB) may exceed available memory ({available_memory:.2f} MB)\n")
    with pyhmmer.easel.SequenceFile(args.faa, digital=True) as seq_file:
        sequences = seq_file.read_block()

    # Run hmmsearch
    results = pyhmmer.hmmer.hmmsearch(hmms, sequences, cpus=args.threads)

    # Open output files conditionally
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
            hmm_name = hmm.name.decode("utf-8", errors="ignore")  # Keep .sr in query name
            hmm_acc = hmm.accession.decode("utf-8", errors="ignore") if hmm.accession else "-"
            
            if allresults_file:
                allresults_file.write(f"HMM: {hmm_name} (accession: {hmm_acc})\n")

            for hit in result:
                hit_name = hit.name.decode("utf-8", errors="ignore")

                # Full sequence stats
                full_evalue = hit.evalue
                full_score = hit.score
                full_bias = hit.bias

                # Write to tblout with query accession
                tblout.write(f"{hit_name:<20} {hmm_name:<20} {hmm_acc:<12} {full_evalue:<9.1e} {full_score:<6.1f} {full_bias:<5.1f}\n")

                # Optional outputs
                if allresults_file:
                    hit_acc = hit.accession.decode("utf-8", errors="ignore") if hit.accession else "-"
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

if __name__ == "__main__":
    main()