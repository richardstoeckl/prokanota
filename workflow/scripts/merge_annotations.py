#!/usr/bin/env python3
"""
Merge all database annotation results with the base feature table.
Dynamically handles any number of databases based on config.
"""

import argparse
import polars as pl
from pathlib import Path
import json


def load_base_table(path: Path) -> pl.DataFrame:
    """Load the base annotation table (e.g., from feature prediction)."""
    return pl.read_csv(path, separator="\t")


def load_db_annotation(path: Path) -> pl.DataFrame:
    """Load a single database annotation result."""
    return pl.read_csv(path, separator="\t")


def merge_annotations(
    base_df: pl.DataFrame,
    db_annotations: list[tuple[str, Path]],
    gene_id_column: str = "gene_id",
) -> pl.DataFrame:
    """
    Merge base table with all database annotations.

    Args:
        base_df: Base annotation dataframe with gene features
        db_annotations: List of (db_name, path) tuples, sorted by order
        gene_id_column: Column name to join on
    """
    result = base_df

    for db_name, db_path in db_annotations:
        if not db_path.exists():
            print(f"Warning: {db_path} does not exist, skipping {db_name}")
            continue

        db_df = load_db_annotation(db_path)

        # Get columns to add (everything except gene_id)
        new_cols = [c for c in db_df.columns if c != gene_id_column]

        # Left join to preserve all genes
        result = result.join(
            db_df,
            on=gene_id_column,
            how="left",
            suffix=f"_{db_name}",  # Prevent collisions
        )

        # Fill nulls with "*" for missing annotations
        result = result.with_columns([
            pl.col(c).cast(pl.Utf8).fill_null("*") for c in new_cols
        ])

    return result


def main():
    parser = argparse.ArgumentParser(description="Merge all annotations into final table")
    parser.add_argument("--base-table", required=True, help="Path to base feature table")
    parser.add_argument("--output", required=True, help="Path to output final annotation TSV")
    parser.add_argument(
        "--db-results",
        required=True,
        help='JSON array of {"name": "db_name", "path": "/path/to/results.tsv", "order": 10}',
    )
    parser.add_argument("--gene-id-column", default="gene_id", help="Column name for gene ID")
    args = parser.parse_args()

    db_results = json.loads(args.db_results)

    # Sort by order
    db_results_sorted = sorted(db_results, key=lambda x: x["order"])

    # Load base table
    base_df = load_base_table(Path(args.base_table))

    # Prepare list of (name, path) tuples
    db_annotations = [(db["name"], Path(db["path"])) for db in db_results_sorted]

    # Merge
    final_df = merge_annotations(base_df, db_annotations, args.gene_id_column)

    # Write output
    final_df.write_csv(args.output, separator="\t")
    print(f"Final annotation table written to {args.output}")
    print(f"Shape: {final_df.shape[0]} rows x {final_df.shape[1]} columns")


if __name__ == "__main__":
    main()
