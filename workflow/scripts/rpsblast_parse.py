#!/usr/bin/env python3
"""
Parse RPS-BLAST tabular output and join with mapping file.
Uses Polars for efficient data processing.
"""

import argparse
import logging
import polars as pl
from pathlib import Path

logger = logging.getLogger(__name__)


def parse_rpsblast_output(output_path: Path) -> pl.DataFrame:
    """
    Parse RPS-BLAST outfmt 6 output.
    Columns: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    """
    if output_path.stat().st_size == 0:
        logger.warning("RPS-BLAST output '%s' is empty; returning no hits", output_path)
        return pl.DataFrame(
            schema={
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
            }
        )

    return pl.read_csv(
        output_path,
        separator="\t",
        has_header=False,
        new_columns=[
            "gene_id", "accession", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "score"
        ],
    )
    # ).with_columns(
    #     pl.col("accession").str.replace("^gnl\\|CDD\\|", "").str.replace("^CDD:", "") # TODO this should not be here, it should be part of the mapping parsing
    # )


def parse_mapping(mapping_path: Path) -> pl.DataFrame:
    """
    Parse mapping file with schema:
    accession<TAB>short_name<TAB>description<TAB>category
    """
    if mapping_path.stat().st_size == 0:
        logger.warning("Mapping file '%s' is empty; returning no annotations", mapping_path)
        return pl.DataFrame(
            schema={
                "accession": pl.Utf8,
                "short_name": pl.Utf8,
                "description": pl.Utf8,
                "category": pl.Utf8,
            }
        )

    return pl.read_csv(
        mapping_path,
        separator="\t",
        has_header=False,
        new_columns=["accession", "short_name", "description", "category"],
        schema_overrides=[pl.Utf8, pl.Utf8, pl.Utf8, pl.Utf8],
    )


def filter_and_deduplicate(df: pl.DataFrame, evalue_cutoff: float = 1e-3) -> pl.DataFrame:
    """
    Apply e-value cutoff, then keep best hit per gene_id.
    Best = highest score, then lowest e-value (ties).
    """
    return (
        df.filter(pl.col("evalue") <= evalue_cutoff)
        .sort(["score", "evalue"], descending=[True, False])
        .unique(subset=["gene_id"], keep="first")
    )


def join_with_mapping(hits_df: pl.DataFrame, mapping_df: pl.DataFrame) -> pl.DataFrame:
    """
    Join hits with mapping on accession.
    """
    return hits_df.join(
        mapping_df,
        left_on="accession",
        right_on="accession",
        how="left",
    )


def format_output(df, db_name, columns_config):
    """Select and rename columns according to the database config."""
    exprs = [pl.col("gene_id")]
    for col_def in columns_config:
        source = col_def["source"]
        col_name = col_def["name"]
        if source in df.columns:
            exprs.append(pl.col(source).alias(col_name))
        else:
            exprs.append(pl.lit(None).alias(col_name))
    return df.select(exprs).fill_null("*")


def main():
    parser = argparse.ArgumentParser(description="Parse RPS-BLAST results and join with mapping")
    parser.add_argument("--rpsblast-output", required=True, help="Path to RPS-BLAST tabular output file")
    parser.add_argument("--mapping", required=True, help="Path to mapping TSV file")
    parser.add_argument("--output", required=True, help="Path to output TSV file")
    parser.add_argument("--db-name", required=True, help="Database name for column prefixes")
    parser.add_argument("--evalue-cutoff", type=float, default=1e-3, help="E-value cutoff")
    parser.add_argument(
        "--columns",
        required=True,
        help="Column config as JSON string",
    )
    args = parser.parse_args()

    import json
    columns_config = json.loads(args.columns)

    # Parse inputs
    hits_df = parse_rpsblast_output(Path(args.rpsblast_output))
    mapping_df = parse_mapping(Path(args.mapping))

    # Process
    filtered_df = filter_and_deduplicate(hits_df, args.evalue_cutoff)
    joined_df = join_with_mapping(filtered_df, mapping_df)
    output_df = format_output(joined_df, args.db_name, columns_config)

    # Write output
    output_df.write_csv(args.output, separator="\t")


if __name__ == "__main__":
    main()
