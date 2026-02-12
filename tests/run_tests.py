#!/usr/bin/env python3
"""Run a minimal Snakemake test run and compare outputs to fixtures."""

from __future__ import annotations

import difflib
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = REPO_ROOT / "workflow" / "Snakefile"
CONFIGFILE = REPO_ROOT / "config" / "test-config.yaml"
INTERIM_DIR = REPO_ROOT / "tests" / "output" / "interim"
RESULTS_DIR = REPO_ROOT / "tests" / "output" / "results"
EXPECTED_DIR = REPO_ROOT / "tests" / "expected" / "pyhmmer"
EXPECTED_FEATURES_DIR = REPO_ROOT / "tests" / "expected" / "features"
EXPECTED_ANNOTATION_DIR = REPO_ROOT / "tests" / "expected" / "annotation"
DB_NAME = "test_pyhmmer"
SAMPLES = [
    "Pyrococcus_furiosus_DSM_3638-GCF_008245085",
    "Pyrococcus_furiosus_DSM_3638-GCF_008245085_protein",
]


def run_snakemake() -> None:
    cmd = [
        "snakemake",
        "--snakefile",
        str(SNAKEFILE),
        "--configfile",
        str(CONFIGFILE),
        "--cores",
        "1",
        "--forceall",
        "--sdm",
        "conda",
    ]
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, cwd=REPO_ROOT)
    if result.returncode != 0:
        sys.exit(result.returncode)


def normalize_tblout(path: Path) -> str:
    if not path.exists():
        raise FileNotFoundError(f"Missing tblout: {path}")
    lines = []
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        lines.append(" ".join(parts))
    return "\n".join(lines).strip()


def normalize_parsed(path: Path) -> str:
    if not path.exists():
        raise FileNotFoundError(f"Missing parsed TSV: {path}")
    lines = [line.rstrip("\n") for line in path.read_text().splitlines() if line.strip()]
    if not lines:
        return ""
    header = lines[0]
    rows = [line.split("\t") for line in lines[1:]]
    if not rows:
        return header.strip()
    header_cols = header.split("\t")
    try:
        gene_idx = header_cols.index("gene_id")
    except ValueError:
        gene_idx = 0
    rows_sorted = sorted(rows, key=lambda r: r[gene_idx])
    normalized_rows = ["\t".join(row) for row in rows_sorted]
    return "\n".join([header] + normalized_rows).strip()


def compare_files(actual: Path, expected: Path, normalizer, label: str) -> bool:
    if not expected.exists():
        print(f"Missing expected fixture for {label}: {expected}")
        return False
    try:
        actual_text = normalizer(actual)
        expected_text = normalizer(expected)
    except FileNotFoundError as exc:
        print(str(exc))
        return False

    if actual_text == expected_text:
        print(f"OK: {label}")
        return True

    diff = difflib.unified_diff(
        expected_text.splitlines(),
        actual_text.splitlines(),
        fromfile=str(expected),
        tofile=str(actual),
        lineterm="",
    )
    print(f"Mismatch in {label}:")
    for line in diff:
        print(line)
    return False


def compare_files_bytes(actual: Path, expected: Path, label: str) -> bool:
    if not expected.exists():
        print(f"Missing expected fixture for {label}: {expected}")
        return False
    if not actual.exists():
        print(f"Missing actual output for {label}: {actual}")
        return False

    actual_bytes = actual.read_bytes()
    expected_bytes = expected.read_bytes()
    if actual_bytes == expected_bytes:
        print(f"OK: {label}")
        return True

    actual_text = actual_bytes.decode("utf-8", errors="replace")
    expected_text = expected_bytes.decode("utf-8", errors="replace")
    diff = difflib.unified_diff(
        expected_text.splitlines(),
        actual_text.splitlines(),
        fromfile=str(expected),
        tofile=str(actual),
        lineterm="",
    )
    print(f"Mismatch in {label}:")
    for line in diff:
        print(line)
    return False


def compare_directory(expected_root: Path, actual_root: Path, label: str) -> bool:
    if not expected_root.exists():
        print(f"Missing expected fixture directory for {label}: {expected_root}")
        return False
    if not actual_root.exists():
        print(f"Missing actual output directory for {label}: {actual_root}")
        return False

    expected_files = {path.relative_to(expected_root) for path in expected_root.rglob("*") if path.is_file()}
    actual_files = {path.relative_to(actual_root) for path in actual_root.rglob("*") if path.is_file()}
    success = True

    missing = sorted(expected_files - actual_files)
    extra = sorted(actual_files - expected_files)
    if missing:
        print(f"Missing expected files for {label}:")
        for path in missing:
            print(f"  {path}")
        success = False
    if extra:
        print(f"Unexpected extra files for {label}:")
        for path in extra:
            print(f"  {path}")
        success = False

    for rel_path in sorted(expected_files & actual_files):
        expected_file = expected_root / rel_path
        actual_file = actual_root / rel_path
        if not compare_files_bytes(actual_file, expected_file, f"{label}:{rel_path}"):
            success = False

    return success


def main() -> int:
    run_snakemake()

    success = True
    for sample_id in SAMPLES:
        actual_tblout = INTERIM_DIR / sample_id / DB_NAME / "hits.tblout"
        expected_tblout = EXPECTED_DIR / sample_id / "hits.tblout"
        if not compare_files(actual_tblout, expected_tblout, normalize_tblout, f"raw:{sample_id}"):
            success = False

        actual_parsed = INTERIM_DIR / sample_id / DB_NAME / "parsed_annotation.tsv"
        expected_parsed = EXPECTED_DIR / sample_id / "parsed_annotation.tsv"
        if not compare_files(actual_parsed, expected_parsed, normalize_parsed, f"parsed:{sample_id}"):
            success = False

        expected_features = EXPECTED_FEATURES_DIR / sample_id
        actual_features = RESULTS_DIR / sample_id / "features"
        if not compare_directory(expected_features, actual_features, f"features:{sample_id}"):
            success = False

        expected_annotation = EXPECTED_ANNOTATION_DIR / sample_id
        actual_annotation = RESULTS_DIR / sample_id / "annotation"
        if not compare_directory(expected_annotation, actual_annotation, f"annotation:{sample_id}"):
            success = False

    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
