#!/usr/bin/env python3
"""Unit tests for mapping file parsing and validation utilities."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = REPO_ROOT / "workflow" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from mapping_utils import parse_mapping_file, validate_mapping_file


class MappingUtilsTests(unittest.TestCase):
    """Coverage for robust mapping parsing behavior."""

    def _write_temp_mapping(self, content: str) -> Path:
        temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(temp_dir.cleanup)
        path = Path(temp_dir.name) / "mapping.tsv"
        path.write_text(content, encoding="utf-8")
        return path

    def _records(self, mapping_df) -> list[dict[str, str]]:
        if hasattr(mapping_df, "to_dicts"):
            return mapping_df.to_dicts()
        return mapping_df.to_dict(orient="records")

    def test_valid_mapping_parses(self) -> None:
        mapping_path = self._write_temp_mapping(
            "ACC001\tGeneA\tDescription A\tcatA\n"
            "ACC002\tGeneB\tDescription B\tcatB\n"
        )

        mapping_df = parse_mapping_file(mapping_path)
        records = self._records(mapping_df)

        self.assertEqual(len(records), 2)
        self.assertEqual(records[0]["accession"], "ACC001")
        self.assertEqual(records[1]["short_name"], "GeneB")

    def test_empty_markers_are_normalized_to_asterisk(self) -> None:
        mapping_path = self._write_temp_mapping(
            "ACC001\t\tNA\tNULL\n"
            "ACC002\t-\tN/A\t*\n"
            "ACC003\t  \tnull\t   \n"
        )

        mapping_df = parse_mapping_file(mapping_path)
        records = self._records(mapping_df)

        for row in records:
            self.assertEqual(row["short_name"], "*")
            self.assertEqual(row["description"], "*")
            self.assertEqual(row["category"], "*")

    def test_wrong_column_count_raises(self) -> None:
        mapping_path = self._write_temp_mapping(
            "ACC001\tGeneA\tDescA\tcatA\n"
            "ACC002\tGeneB\tDescB\n"
        )

        with self.assertRaises(ValueError) as ctx:
            parse_mapping_file(mapping_path)

        self.assertIn("Expected 4 columns, found 3", str(ctx.exception))

    def test_duplicate_accession_raises(self) -> None:
        mapping_path = self._write_temp_mapping(
            "ACC001\tGeneA\tDescA\tcatA\n"
            "ACC001\tGeneB\tDescB\tcatB\n"
        )

        with self.assertRaises(ValueError) as ctx:
            parse_mapping_file(mapping_path)

        self.assertIn("Duplicate accession 'ACC001'", str(ctx.exception))

    def test_missing_accession_raises(self) -> None:
        mapping_path = self._write_temp_mapping(
            "\tGeneA\tDescA\tcatA\n"
            "ACC002\tGeneB\tDescB\tcatB\n"
        )

        with self.assertRaises(ValueError) as ctx:
            parse_mapping_file(mapping_path)

        self.assertIn("Column 1 (accession) must not be empty", str(ctx.exception))

    def test_special_characters_in_description_are_preserved(self) -> None:
        mapping_path = self._write_temp_mapping(
            "ACC001\tGeneA\tcontains 'single' and \"double\" quotes Ω\tcatA\n"
        )

        mapping_df = parse_mapping_file(mapping_path)
        records = self._records(mapping_df)

        self.assertEqual(records[0]["description"], "contains 'single' and \"double\" quotes Ω")

    def test_wrapped_double_quotes_are_stripped(self) -> None:
        mapping_path = self._write_temp_mapping(
            '"CDD:XXXXX"\t"gene_1"\t"description with spaces"\t"category"\n'
        )

        mapping_df = parse_mapping_file(mapping_path)
        records = self._records(mapping_df)

        self.assertEqual(records[0]["accession"], "CDD:XXXXX")
        self.assertEqual(records[0]["short_name"], "gene_1")
        self.assertEqual(records[0]["description"], "description with spaces")
        self.assertEqual(records[0]["category"], "category")

    def test_validate_mapping_file_returns_stats(self) -> None:
        mapping_path = self._write_temp_mapping(
            "ACC001\tGeneA\tDescA\tcatA\n"
            "ACC002\tGeneB\tDescB\tcatB\n"
            "ACC003\tGeneC\tDescC\tcatC\n"
        )

        stats = validate_mapping_file(mapping_path)

        self.assertEqual(stats["row_count"], 3)
        self.assertEqual(stats["sample_accessions"], ["ACC001", "ACC002", "ACC003"])


if __name__ == "__main__":
    unittest.main()
