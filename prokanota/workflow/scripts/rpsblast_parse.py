#!/usr/bin/env python3
"""Read rpsblast hits and pass them to shared annotation processing."""

from pathlib import Path

import polars as pl
from parse_utils import filter_and_deduplicate as filter_and_deduplicate
from parse_utils import run_parser


def parse_rpsblast_output(output_path: Path) -> pl.DataFrame:
    """
    Parse RPS-BLAST outfmt 6 output.
    Columns: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen
    """
    schema = {
        "gene_id": pl.Utf8,
        "accession": pl.Utf8,
        "pident": pl.Float64,
        "length": pl.Int64,
        "mismatch": pl.Int64,
        "gapopen": pl.Int64,
        "qstart": pl.Int64,
        "qend": pl.Int64,
        "sstart": pl.Int64,
        "send": pl.Int64,
        "evalue": pl.Float64,
        "score": pl.Float64,
        "subject_length": pl.Int64,
    }
    if output_path.stat().st_size == 0:
        return pl.DataFrame(schema=schema)

    return pl.read_csv(
        output_path,
        separator="\t",
        has_header=False,
        schema=schema,
    )


def main():
    run_parser("RPS-BLAST", "--rpsblast-output", parse_rpsblast_output)


if __name__ == "__main__":
    main()
