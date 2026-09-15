#!/usr/bin/env python3
"""Read pyhmmer hits and pass them to shared annotation processing."""

import logging
from pathlib import Path

import polars as pl
from parse_utils import run_parser

logger = logging.getLogger(__name__)


def parse_pyhmmer_tblout(tblout_path: Path) -> pl.DataFrame:
    """
    Parse pyhmmer tblout format:
    # target name        query name           accession    E-value  score  bias  model coverage
    """
    # do I really need to check the file size and the "if not rows" later again or does one check suffice?
    # -> The first check catches truly empty files upfront (fast). The second catches files with only comments or whitespace (which would parse to zero rows).
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
                "reference_coverage": pl.Float64,
            }
        )
    rows = []
    with open(tblout_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 7:
                rows.append(
                    {
                        "gene_id": parts[0],
                        "query_name": parts[1],
                        "accession": parts[2],
                        "evalue": float(parts[3]),
                        "score": float(parts[4]),
                        "bias": float(parts[5]),
                        "reference_coverage": float(parts[6]),
                    }
                )
    if not rows:
        logger.warning(
            "pyhmmer tblout '%s' contains no parsable rows; returning no hits",
            tblout_path,
        )
        return pl.DataFrame(
            schema={
                "gene_id": pl.Utf8,
                "query_name": pl.Utf8,
                "accession": pl.Utf8,
                "evalue": pl.Float64,
                "score": pl.Float64,
                "bias": pl.Float64,
                "reference_coverage": pl.Float64,
            }
        )
    return pl.DataFrame(rows)


def filter_and_deduplicate(
    df: pl.DataFrame, evalue_cutoff: float = 1e-3
) -> pl.DataFrame:
    """
    Apply e-value cutoff, then keep best hit per gene_id.
    Best = highest score, lowest e-value, highest model coverage, then accession
    (or query name when the accession is missing).
    """
    return (
        df.filter(pl.col("evalue") <= evalue_cutoff)
        .with_columns(
            pl.when(pl.col("accession") == "-")
            .then(pl.col("query_name"))
            .otherwise(pl.col("accession"))
            .alias("identifier")
        )
        .sort(
            ["score", "evalue", "reference_coverage", "identifier"],
            descending=[True, False, True, False],
        )
        .unique(subset=["gene_id"], keep="first")
    )


def format_output(df, db_name, columns_config):
    """Select and rename columns according to the database config."""
    exprs = [pl.col("gene_id")]
    for col_def in columns_config:
        source = col_def["source"]
        col_name = col_def["name"]
        if source in df.columns:
            expression = pl.col(source)
            if source in ("score", "bias"):
                expression = expression.round(1)
            elif source == "evalue":
                expression = expression.map_elements(
                    lambda value: float(f"{value:.1e}"),
                    return_dtype=pl.Float64,
                )
            exprs.append(expression.alias(col_name))
        else:
            exprs.append(pl.lit(None).alias(col_name))
    return df.select(exprs).fill_null("*")


def main():
    run_parser(
        "pyhmmer",
        "--tblout",
        parse_pyhmmer_tblout,
        filter_and_deduplicate,
        format_output,
    )


if __name__ == "__main__":
    main()
