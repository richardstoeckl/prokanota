#!/usr/bin/env python3
"""
Copyright Richard Stöckl 2025.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE or copy at 
https://www.boost.org/LICENSE_1_0.txt)

Steps:
  1. Extract the tar.gz archive containing .sr files, stripping the top-level directory.
  2. Build HMMs concurrently from each .sr file using a thread-local pyhmmer.plan7.Builder.
  3. Press the list of HMMs into a combined HMM database using pyhmmer.hmmer.hmmpress.

Usage:
    python prepareArCOGdb.py --tar <path_to_tar.gz> --output-dir <output_directory> [--workers 4]
"""

import argparse
import os
import tarfile
from glob import glob
from concurrent.futures import ThreadPoolExecutor, as_completed
import traceback

import pyhmmer.easel
import pyhmmer.plan7
import pyhmmer.hmmer

def extract_tar_strip(tar_path, extract_dir):
    """Extract the tar.gz archive into extract_dir, stripping the top-level directory."""
    os.makedirs(extract_dir, exist_ok=True)
    with tarfile.open(tar_path, "r:gz") as tar:
        for member in tar.getmembers():
            # Strip the top-level directory
            member.name = os.path.basename(member.name)
            if member.isfile():
                tar.extract(member, path=extract_dir)
    print(f"Extracted files to {extract_dir}")

def process_sr_file(sr_file, alphabet):
    sr_name = os.path.basename(sr_file).replace(".sr", "")
    print(f"Processing {sr_name}")
    try:
        with pyhmmer.easel.MSAFile(sr_file, digital=True, alphabet=alphabet) as msa_file:
            msa = msa_file.read()
            if msa is None or len(list(msa.sequences)) == 0:
                raise ValueError(f"No alignment data in {sr_name}")
    except Exception as e:
        print(f"Error reading MSA from {sr_name}: {e}")
        raise

    # Convert sequences to a list to allow slicing
    sequences = list(msa.sequences)

    # Set the MSA name if missing
    msa.name = sr_name.encode("utf-8")
    if not msa.name:
        raise ValueError(f"Empty name for {sr_name}")
    print(f"MSA name set to: {msa.name}")

    # Set the MSA accession if missing
    if msa.accession is None:
        msa.accession = sr_name.encode("utf-8")
        if not msa.accession:
            raise ValueError(f"Empty accession for {sr_name}")

    # Check if all sequences are identical by comparing their string representations
    first_seq_str = str(sequences[0])
    if all(str(seq) == first_seq_str for seq in sequences[1:]):
        print(f"Warning: All sequences in {sr_name} are identical. Skipping.")
        return None

    print(f"MSA {sr_name}: {len(sequences)} sequences, alignment length {len(sequences[0])}")

    # Create thread-local Builder and Background instances
    builder = pyhmmer.plan7.Builder(alphabet=alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    
    try:
        hmm, _, _ = builder.build_msa(msa, background)
    except Exception as e:
        print(f"Error building HMM for {sr_name}: {e}")
        raise
    return hmm

def build_hmms_concurrently(extract_dir, alphabet, max_workers=4):
    """Build HMMs concurrently from all .sr files in extract_dir."""
    sr_files = glob(os.path.join(extract_dir, "*.sr"))
    if not sr_files:
        raise ValueError(f"No .sr files found in {extract_dir}")

    hmms = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_sr_file, sr, alphabet): sr for sr in sr_files}
        for future in as_completed(futures):
            sr_path = futures[future]
            try:
                hmm = future.result()
                if hmm is not None:
                    hmms.append(hmm)
            except Exception as e:
                print(f"Error processing {sr_path}: {e}")
                traceback.print_exc()
    if not hmms:
        raise ValueError("No HMMs were successfully built")
    return hmms

def press_hmm_db(hmms, output):
    """
    Press the list of HMMs into a database.
    This will create the files:
      {output}.h3p, {output}.h3m, {output}.h3f, {output}.h3i
    """
    pyhmmer.hmmer.hmmpress(hmms, output)
    print(f"Pressed HMM database created for {output}")

def main():
    parser = argparse.ArgumentParser(
        description="Prepare arCOG HMM database using pyhmmer with concurrent HMM building")
    parser.add_argument("--tar", required=True, help="Path to the downloaded arCOG tar.gz file")
    parser.add_argument("--output-dir", required=True, help="Directory to store the prepared DB")
    parser.add_argument("--workers", type=int, default=4, help="Number of worker threads")
    args = parser.parse_args()

    extract_dir = os.path.join(args.output_dir, "arCOG_raw")
    # The output path here will serve as the base for the pressed database files
    output_base = os.path.join(args.output_dir, "arCOGs_combined")
    alphabet = pyhmmer.easel.Alphabet.amino()

    # Step 1: Extract the tar.gz archive
    extract_tar_strip(args.tar, extract_dir)

    # Step 2: Build HMMs concurrently (each thread instantiates its own Builder)
    hmms = build_hmms_concurrently(extract_dir, alphabet, max_workers=args.workers)

    # Step 3: Press the combined HMM database to create index files
    press_hmm_db(hmms, output_base)

if __name__ == "__main__":
    main()
