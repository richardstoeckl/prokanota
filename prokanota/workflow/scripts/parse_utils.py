"""Shared annotation processing after tool-specific search-result reading."""

import argparse
import json
from datetime import datetime
from pathlib import Path

import polars as pl
from logging_utils import (
    build_part,
    log_file_paths,
    log_parameters,
    log_prokanota_version,
    log_statistics,
    setup_logger,
)
from mapping_utils import parse_mapping_file


def filter_and_deduplicate(
    df: pl.DataFrame, evalue_cutoff: float = 1e-3
) -> pl.DataFrame:
    """
    Apply e-value cutoff, then keep best hit per gene_id.
    Best = highest score, lowest e-value, highest subject coverage, then accession.
    """
    return (
        df.filter(pl.col("evalue") <= evalue_cutoff)
        .with_columns(
            ((pl.col("send") - pl.col("sstart")).abs() + 1)
            .truediv(pl.col("subject_length"))
            .alias("reference_coverage")
        )
        .sort(
            ["score", "evalue", "reference_coverage", "accession"],
            descending=[True, False, True, False],
        )
        .unique(subset=["gene_id"], keep="first")
    )


def join_hits_with_mapping(
    hits_df: pl.DataFrame, mapping_df: pl.DataFrame, join_on: str = "accession"
) -> pl.DataFrame:
    """
    Join hits with mapping on specified column.

    Args:
        hits_df: DataFrame with gene_id, accession, and other hit columns
        mapping_df: DataFrame with accession (key), short_name, description, category
        join_on: Column name to join on. Valid values: "accession" or "query_name".
                The mapping file must have the same values in its accession column as the selected join_on column.
                No automatic normalization is performed; keys must match exactly.

    Returns:
        Left-joined DataFrame with all columns from hits and mapping annotations.
        Unmatched hits will have NA values for mapping columns.
    """
    # Validate join_on parameter
    if join_on not in ("accession", "query_name"):
        raise ValueError(
            f"join_on must be 'accession' or 'query_name', got '{join_on}'"
        )

    if join_on not in hits_df.columns:
        raise ValueError(
            f"mapping_key '{join_on}' is not present in search results columns: {hits_df.columns}"
        )

    return hits_df.join(
        mapping_df,
        left_on=join_on,
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


def run_parser(
    tool_name,
    input_option,
    read_hits,
    select_hits=filter_and_deduplicate,
    format_hits=format_output,
):
    """Run shared filtering, mapping and output logic with explicit tool policies."""
    parser = argparse.ArgumentParser(
        description=f"Parse {tool_name} results and join with mapping"
    )
    parser.add_argument(
        input_option,
        dest="search_output",
        required=True,
        help=f"Path to {tool_name} search output",
    )
    parser.add_argument("--mapping", required=True, help="Path to mapping TSV file")
    parser.add_argument("--output", required=True, help="Path to output TSV file")
    parser.add_argument(
        "--db-name", required=True, help="Database name for column prefixes"
    )
    parser.add_argument(
        "--evalue-cutoff", type=float, default=1e-3, help="E-value cutoff"
    )
    parser.add_argument(
        "--mapping-key",
        type=str,
        default="accession",
        choices=["accession", "query_name"],
        help="Column to join on with mapping file (default: accession)",
    )
    parser.add_argument(
        "--columns",
        required=True,
        help="Column config as JSON string",
    )
    args = parser.parse_args()

    logger = setup_logger(build_part("PARSE", args.db_name))
    log_prokanota_version(logger)
    start_time = datetime.now()
    logger.info("=" * 60)
    logger.info(f"{tool_name} Parse Started")
    logger.info("=" * 60)

    # Log input files
    log_file_paths(
        logger,
        search_output=args.search_output,
        mapping_file=args.mapping,
        output_file=args.output,
    )

    # Log parameters
    log_parameters(
        logger,
        database=args.db_name,
        evalue_cutoff=args.evalue_cutoff,
    )

    columns_config = json.loads(args.columns)

    # Parse inputs
    logger.info(f"Parsing {tool_name} search results...")
    hits_df = read_hits(Path(args.search_output))
    hits_count = len(hits_df)
    if hits_count == 0:
        logger.warning("Search output %s contains no hits", args.search_output)
    logger.info(f"Search results rows read: {hits_count}")

    logger.info("Parsing mapping file...")
    mapping_df = parse_mapping_file(Path(args.mapping))
    if mapping_df.height == 0:
        logger.warning(
            "Mapping file %s is empty; returning no annotations", args.mapping
        )
    mapping_count = len(mapping_df)
    logger.info(f"Mapping file rows read: {mapping_count}")

    # Process
    logger.info("Filtering by e-value and selecting unique top hits...")
    passing_count = hits_df.filter(pl.col("evalue") <= args.evalue_cutoff).height
    filtered_df = select_hits(hits_df, args.evalue_cutoff)
    filtered_count = len(filtered_df)
    filtered_out = hits_count - passing_count
    log_statistics(
        logger,
        input_hits=hits_count,
        rejected_by_evalue=filtered_out,
        passing_evalue=passing_count,
        alternative_hits_removed=passing_count - filtered_count,
        unique_top_hits=filtered_count,
    )

    # Check for excessive filtering
    if hits_count > 0 and filtered_out > (hits_count * 0.95):
        logger.warning(
            f"{100 * filtered_out / hits_count:.1f}% of hits were filtered by e-value cutoff. Consider checking parameters or cutoff threshold."
        )

    logger.info(f"Joining with mapping file using mapping_key='{args.mapping_key}'...")
    joined_df = join_hits_with_mapping(
        filtered_df, mapping_df, join_on=args.mapping_key
    )
    matched_count = joined_df.filter(pl.col("short_name").is_not_null()).shape[0]
    unmatched_count = filtered_count - matched_count
    log_statistics(
        logger,
        matched_with_mapping=matched_count,
        unmatched_with_mapping=unmatched_count,
        mapping_key=args.mapping_key,
    )

    # Check mapping match results
    if filtered_count > 0:
        if matched_count == 0:
            example_hit = (
                filtered_df[args.mapping_key][0] if filtered_count > 0 else "<none>"
            )
            example_map = mapping_df["accession"][0] if mapping_count > 0 else "<none>"
            logger.warning(
                f"NO HITS matched the mapping file! Check if accession format is correct in search results vs mapping file. "
                f"Here, hits with an accession format like '{example_hit}' were unable to be joined with {args.mapping_key}s from the mapping file like '{example_map}'."
            )
        elif matched_count < (filtered_count * 0.5):
            logger.warning(
                f"Low mapping match rate - only {matched_count}/{filtered_count} hits ({100 * matched_count / filtered_count:.1f}%) matched. Check accession format consistency."
            )

    logger.info("Formatting output...")
    for col_def in columns_config:
        if col_def["source"] not in joined_df.columns:
            raise ValueError(
                f"Database '{args.db_name}' ({tool_name}): output column "
                f"'{col_def['name']}' requests unavailable source '{col_def['source']}'. "
                f"Available sources: {', '.join(joined_df.columns)}"
            )
    output_df = format_hits(joined_df, args.db_name, columns_config).sort("gene_id")

    # Write output
    logger.info(f"Writing output to {args.output}...")
    output_df.write_csv(args.output, separator="\t")
    output_count = len(output_df)
    log_statistics(logger, output_rows=output_count)

    # Check for empty output if we had input
    if output_count == 0 and hits_count > 0:
        logger.warning(
            f"All {hits_count} input hits were lost during processing. Check data consistency or filtering parameters."
        )

    # Log execution time
    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    logger.info(f"Execution time: {duration:.2f} seconds")
    logger.info("=" * 60)
    logger.info(f"{tool_name} Parse Completed Successfully")
    logger.info("=" * 60)
