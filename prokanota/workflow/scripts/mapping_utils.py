"""
Shared mapping-file parsing and validation utilities.
"""

from pathlib import Path
from typing import Any

import pandas as pd

MAPPING_COLUMNS = ["accession", "short_name", "description", "category"]
CANONICAL_EMPTY_VALUE = "*"
EMPTY_MARKERS = {"", "*", "-", "na", "n/a", "null"}


def _to_parser_dataframe(df: pd.DataFrame) -> Any:
    """
    Return a DataFrame compatible with parser scripts.

    If polars is available (parse-script environments), return a Polars DataFrame.
    Otherwise (e.g., Snakefile validation in base Snakemake env), return pandas DataFrame.
    """
    try:
        import polars as pl

        return pl.DataFrame({col: df[col].tolist() for col in MAPPING_COLUMNS})
    except ImportError:
        return df


def _empty_mapping_df() -> Any:
    return _to_parser_dataframe(pd.DataFrame(columns=MAPPING_COLUMNS))


def _normalize_mapping_value(value: str | None) -> str:
    if value is None:
        return CANONICAL_EMPTY_VALUE

    stripped = _strip_wrapping_quotes(value)
    if stripped.lower() in EMPTY_MARKERS:
        return CANONICAL_EMPTY_VALUE
    return stripped


def _strip_wrapping_quotes(value: str) -> str:
    stripped = value.strip()
    if len(stripped) >= 2 and stripped.startswith('"') and stripped.endswith('"'):
        return stripped[1:-1]
    return stripped


def _is_missing_key(value: str | None) -> bool:
    if value is None:
        return True

    stripped = _strip_wrapping_quotes(value)
    return stripped.lower() in EMPTY_MARKERS


def parse_mapping_file(mapping_path: str | Path) -> Any:
    """
    Parse mapping TSV file with schema:
    accession<TAB>short_name<TAB>description<TAB>category

    Rules:
    - No header
    - Exactly 4 tab-separated columns per non-empty line
    - Empty lines are ignored
    - Key column (accession) must not be empty or an empty marker
    - Duplicate accessions are rejected
    - Non-key empty markers are normalized to '*'

    Raises:
        FileNotFoundError: If file does not exist.
        ValueError: If file content is malformed.
    """
    mapping_path = Path(mapping_path)
    if not mapping_path.exists():
        raise FileNotFoundError(f"Mapping file not found: {mapping_path}")

    malformed_rows: list[str] = []
    seen_keys: dict[str, int] = {}
    records: list[dict[str, str]] = []

    with open(mapping_path, "r", encoding="utf-8", errors="replace") as handle:
        for line_num, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\n\r")

            if not line.strip():
                continue

            columns = line.split("\t")
            if len(columns) != 4:
                malformed_rows.append(
                    f"Line {line_num}: Expected 4 columns, found {len(columns)}"
                )
                continue

            accession_raw, short_name_raw, description_raw, category_raw = columns
            accession = _strip_wrapping_quotes(accession_raw)

            if _is_missing_key(accession):
                malformed_rows.append(
                    f"Line {line_num}: Column 1 (accession) must not be empty or a missing marker"
                )
                continue

            if accession in seen_keys:
                malformed_rows.append(
                    f"Line {line_num}: Duplicate accession '{accession}' (first seen at line {seen_keys[accession]})"
                )
                continue

            seen_keys[accession] = line_num
            records.append(
                {
                    "accession": accession,
                    "short_name": _normalize_mapping_value(short_name_raw),
                    "description": _normalize_mapping_value(description_raw),
                    "category": _normalize_mapping_value(category_raw),
                }
            )

    if malformed_rows:
        error_msg = f"Mapping file {mapping_path} has invalid structure:\n" + "\n".join(
            malformed_rows[:10]
        )
        if len(malformed_rows) > 10:
            error_msg += f"\n... and {len(malformed_rows) - 10} more errors"
        raise ValueError(error_msg)

    if not records:
        return _empty_mapping_df()

    mapping_df = pd.DataFrame.from_records(records, columns=MAPPING_COLUMNS)
    return _to_parser_dataframe(mapping_df)


def validate_mapping_file(mapping_path: str | Path) -> dict[str, Any]:
    """
    Validate a mapping file using the shared parser and return summary stats.
    """
    mapping_df = parse_mapping_file(mapping_path)

    if hasattr(mapping_df, "height"):
        row_count = mapping_df.height
        sample_accessions = mapping_df["accession"].head(5).to_list()
    else:
        row_count = len(mapping_df)
        sample_accessions = mapping_df["accession"].head(5).tolist()

    return {
        "row_count": row_count,
        "sample_accessions": sample_accessions,
    }
