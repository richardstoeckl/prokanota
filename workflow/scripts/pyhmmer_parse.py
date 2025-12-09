#!/usr/bin/env python3
"""
Parse pyhmmer tblout results and join with mapping file.
Replicates the shell parsing logic from parseSearchResults.smk using Polars.
"""

import argparse
import logging
import polars as pl
from pathlib import Path

logger = logging.getLogger(__name__)


def parse_pyhmmer_tblout(tblout_path: Path) -> pl.DataFrame:
    """
    Parse pyhmmer tblout format:
    # target name        query name           accession    E-value  score  bias
    """
    # TODO do I really need to check the file size and the "if not rows" later again or does one check suffice?
    if tblout_path.stat().st_size == 0:
        logger.warning("pyhmmer tblout '%s' is empty; returning no hits", tblout_path)
        return pl.DataFrame(
            schema={
                "gene_id": pl.Utf8,
                "query_name": pl.Utf8,
                "accession": pl.Utf8,
                "evalue": pl.Float64,
                "score": pl.Float64,
                "bias": pl.Float64,
            }
        )
    rows = []
    with open(tblout_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 6:
                rows.append({
                    "gene_id": parts[0],
                    "query_name": parts[1],
                    "accession": parts[2],
                    "evalue": float(parts[3]),
                    "score": float(parts[4]),
                    "bias": float(parts[5]),
                })
    if not rows:
        logger.warning("pyhmmer tblout '%s' contains no parsable rows; returning no hits", tblout_path)
        return pl.DataFrame(
            schema={
                "gene_id": pl.Utf8,
                "query_name": pl.Utf8,
                "accession": pl.Utf8,
                "evalue": pl.Float64,
                "score": pl.Float64,
                "bias": pl.Float64,
            }
        )
    return pl.DataFrame(rows)


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
    For pyhmmer results, the 'query_name' or 'accession' column matches the mapping accession.
    """
    # Strip any version suffix (e.g., "arCOG00001.1" -> "arCOG00001")
    # this is outdated
    # TODO: remove if not needed
    # hits_with_clean_acc = hits_df.with_columns(
    #     pl.col("query_name").str.split(".").list.first().alias("query_name_clean")
    # )

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
    parser = argparse.ArgumentParser(description="Parse pyhmmer results and join with mapping")
    parser.add_argument("--tblout", required=True, help="Path to pyhmmer tblout file")
    parser.add_argument("--mapping", required=True, help="Path to mapping TSV file")
    parser.add_argument("--output", required=True, help="Path to output TSV file")
    parser.add_argument("--db-name", required=True, help="Database name for column prefixes")
    parser.add_argument("--evalue-cutoff", type=float, default=1e-3, help="E-value cutoff")
    parser.add_argument(
        "--columns",
        required=True,
        help="Column config as JSON string, e.g., '[{\"name\":\"db_hit\",\"source\":\"query_name\"}]'",
    )
    args = parser.parse_args()

    import json
    columns_config = json.loads(args.columns)

    # Parse inputs
    hits_df = parse_pyhmmer_tblout(Path(args.tblout))
    mapping_df = parse_mapping(Path(args.mapping))

    # Process
    filtered_df = filter_and_deduplicate(hits_df, args.evalue_cutoff)
    joined_df = join_with_mapping(filtered_df, mapping_df)
    output_df = format_output(joined_df, args.db_name, columns_config)

    # Write output
    output_df.write_csv(args.output, separator="\t")


if __name__ == "__main__":
    main()
