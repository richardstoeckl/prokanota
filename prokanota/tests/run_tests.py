#!/usr/bin/env python3
"""Run a minimal Snakemake test run and compare outputs to fixtures."""

from __future__ import annotations

import argparse
import difflib
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
PYTEST_TESTS = PACKAGE_ROOT / "tests" / "test_published_biological_truths.py"
SNAKEFILE = PACKAGE_ROOT / "workflow" / "Snakefile"
CONFIGFILE = PACKAGE_ROOT / "tests" / "test-config.yaml"
DATABASES_FULL = PACKAGE_ROOT / "tests" / "test-databases.yaml"
DATABASES_PYHMMER = PACKAGE_ROOT / "tests" / "test-pyhmmer.yaml"
INTERIM_DIR = PACKAGE_ROOT / "tests" / "output" / "interim"
RESULTS_DIR = PACKAGE_ROOT / "tests" / "output" / "results"
OUTPUT_DIR = PACKAGE_ROOT / "tests" / "output"
EXPECTED_BASE_DIR = PACKAGE_ROOT / "tests" / "expected"
SAMPLES = [
    "Pyrococcus_furiosus_DSM_3638-GCF_008245085",
    "Pyrococcus_furiosus_DSM_3638-GCF_008245085_protein",
]

genome_sample = SAMPLES[0]
crispr_path = RESULTS_DIR / genome_sample / "features" / f"{genome_sample}_crispr.tsv"
rna_path = RESULTS_DIR / genome_sample / "features" / f"{genome_sample}_rna.tsv"
gbk_path = RESULTS_DIR / genome_sample / "features" / f"{genome_sample}.gbk"
gff_path = RESULTS_DIR / genome_sample / "features" / f"{genome_sample}.gff"
tsv_path = RESULTS_DIR / genome_sample / "features" / f"{genome_sample}.tsv"
annotation_path = (
    RESULTS_DIR / genome_sample / "annotation" / f"{genome_sample}_finalAnnotation.tsv"
)
fna_path = RESULTS_DIR / genome_sample / "features" / f"{genome_sample}.fna"
faa_path = RESULTS_DIR / genome_sample / "features" / f"{genome_sample}.faa"


# ANSI color codes for terminal output
class Colors:
    GREEN = "\033[92m"
    RED = "\033[91m"
    YELLOW = "\033[93m"
    BLUE = "\033[94m"
    BOLD = "\033[1m"
    END = "\033[0m"


# To compare the raw outputs of the database searches, we need the corresponding output file.
# Pyhmmer uses the ".tblout" format, all other dbs that are so far implemented just a .tsv format
# so this function swtiches depending on the db.
def get_hits_filename(db_name: str) -> str:
    if db_name == "test_pyhmmer":
        return "hits.tblout"
    return "hits.tsv"


# Raw search outputs from different database tools require specific normalization
# techniques to enable direct, line-by-line file comparisons. This function returns
# the appropriate normalizer function (e.g., normalize_tblout for Pyhmmer) or
# None if no special normalization is needed (e.g., standard TSV tables).
def get_hits_normalizer(db_name: str):
    if db_name == "test_pyhmmer":
        return normalize_tblout
    return None


# Safely loads a YAML configuration file from disk.
# Validates that the file exists, can be parsed by PyYAML, and that the root of
# the YAML document is a dictionary (mapping), raising errors otherwise.
def load_yaml_file(path: Path, label: str) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Missing {label} file: {path}")
    data = yaml.safe_load(path.read_text())
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ValueError(f"Invalid {label} file {path}: root must be a mapping")
    return data


# Loads the database definitions from the registry YAML file.
# Verifies that the file contains a "databases" entry and that its value is a list
# of database configuration mappings.
def load_databases(path: Path) -> list[dict[str, object]]:
    data = load_yaml_file(path, "databases")
    databases = data.get("databases")
    if not isinstance(databases, list):
        raise ValueError(f"Invalid databases file {path}: missing databases list")
    return databases


# Filters a list of database configurations, returning only those that are currently enabled.
# Databases are considered enabled by default unless they have an explicit `enabled: false` key.
def get_enabled_databases(
    databases: list[dict[str, object]],
) -> list[dict[str, object]]:
    enabled = []
    for db in databases:
        if not isinstance(db, dict):
            continue
        if not db.get("enabled", True):
            continue
        enabled.append(db)
    return enabled


# Maps each database name to its defined annotation column names.
# These columns contain specific search/hit metadata (like E-values or scores)
# that are appended to the final annotation TSV file by the workflow.
def get_database_columns(databases: list[dict[str, object]]) -> dict[str, list[str]]:
    columns = {}
    for db in databases:
        name = db.get("name") if isinstance(db, dict) else None
        if not isinstance(name, str):
            continue
        db_columns = []
        for col in db.get("columns", []) if isinstance(db, dict) else []:
            if isinstance(col, dict):
                col_name = col.get("name")
                if isinstance(col_name, str):
                    db_columns.append(col_name)
        columns[name] = db_columns
    return columns


# Collects the set of all annotation columns defined across all databases in the system.
# This set is used to dynamically identify and filter database-related columns
# in the final annotation table during test verification across different execution modes.
def get_all_db_columns(databases: list[dict[str, object]]) -> set[str]:
    all_columns = set()
    for db in databases:
        for col in db.get("columns", []) if isinstance(db, dict) else []:
            if isinstance(col, dict):
                col_name = col.get("name")
                if isinstance(col_name, str):
                    all_columns.add(col_name)
    return all_columns


# Creates a temporary YAML configuration file tailored to the selected test run.
# It reads the base template config, overrides the "databases" setting to point to
# the database YAML for the current mode, and writes the resulting config to the
# log/output directory for Snakemake execution.
def build_snakemake_config(
    config_path: Path, databases_path: Path, mode_label: str
) -> Path:
    config_data = load_yaml_file(config_path, "config")
    global_cfg = config_data.get("global")
    if not isinstance(global_cfg, dict):
        raise ValueError(f"Invalid config file {config_path}: missing global section")
    global_cfg["databases"] = str(databases_path)

    logs_dir = PACKAGE_ROOT / "tests" / "output" / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    target_path = logs_dir / f"test-config-{mode_label}.yaml"
    with target_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(config_data, handle, sort_keys=False)
    return target_path


# Configures and parses the command-line argument parser for the test suite runner.
# Allows selecting the run mode (minimal, pyhmmer, or full) and verbosity settings.
def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run prokanota tests")
    parser.add_argument(
        "--mode",
        choices=("minimal", "pyhmmer", "full"),
        default="pyhmmer",
        help="Select which test mode to run",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show full tool output in the console",
    )
    parser.add_argument(
        "--very-verbose",
        action="store_true",
        help="Show tool output and verbose Snakemake execution details",
    )
    return parser.parse_args(argv)


# Runs Pytest unit tests under test_published_biological_truths.py.
# These unit tests verify biological invariants and pipeline constraints on known reference datasets.
# Returns True if all unit tests passed, False otherwise.
def run_pytest() -> bool:
    cmd = [sys.executable, "-m", "pytest", "-q", str(PYTEST_TESTS)]
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, cwd=PACKAGE_ROOT)
    return result.returncode == 0


# Invokes Snakemake via a subprocess to run the actual Prokanota workflow.
# It forces a complete rerun of all rules (`--forceall`) using the specified configuration
# and isolated Conda environments (`--sdm conda`). Verbosities and print shells are configured
# dynamically based on the passed options.
def run_snakemake(configfile: Path, verbose: bool, very_verbose: bool) -> None:
    if very_verbose:
        verbose = True
    cmd = [
        "snakemake",
        "--snakefile",
        str(SNAKEFILE),
        "--configfile",
        str(configfile),
        "--notemp",
        "--cores",
        "1",
        "--forceall",
        "--sdm",
        "conda",
    ]
    if verbose:
        cmd.extend(["--config", "prokanota_verbose=true"])
    if very_verbose:
        cmd.extend(["--printshellcmds", "--show-failed-logs"])
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, cwd=PACKAGE_ROOT)
    if result.returncode != 0:
        sys.exit(result.returncode)


def parse_fasta(path: Path) -> dict[str, str]:
    """Parse a FASTA file into a dictionary of header -> sequence."""
    if not path.exists():
        raise FileNotFoundError(f"Missing FASTA: {path}")
    records = {}
    current_id = None
    current_seq = []
    # Read and process line by line to handle multi-line fasta entries
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            # If we've already parsed a sequence, save it before starting the next
            if current_id is not None:
                records[current_id] = "".join(current_seq)
            # Extract the identifier (the first whitespace-separated token of the header)
            current_id = line[1:].split()[0] if line[1:].strip() else ""
            current_seq = []
        elif line.strip():
            # Append non-empty sequence chunks
            current_seq.append(line.strip())
    # Save the final sequence entry
    if current_id is not None:
        records[current_id] = "".join(current_seq)
    return records


# A robust utility that parses, filters, sorts, and re-serializes TSV and table files.
# This makes output files order-invariant and allows direct comparison of tables even
# when columns differ across run modes or rows are written out of order.
def normalize_tsv(
    path: Path,
    sort_key_col: str | int = "gene_id",
    column_filter=None,
    delimiter: str = "\t",
    line_parser=None,
    has_header: bool = True,
) -> str:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

    # Strip comments (lines starting with '#') and leading/trailing whitespaces
    raw_lines = path.read_text().splitlines()
    data_lines = [
        line.strip() for line in raw_lines if line.strip() and not line.startswith("#")
    ]
    if not data_lines:
        return ""

    # Split fields using either a custom parser or standard delimiter
    if line_parser is not None:
        rows = [line_parser(line) for line in data_lines]
    else:
        rows = [line.split(delimiter) for line in data_lines]

    if has_header:
        header_cols = rows[0]
        data_rows = rows[1:]
    else:
        header_cols = []
        data_rows = rows

    # Filter columns dynamically if a column filter callback is provided
    if column_filter is not None and header_cols:
        keep_indexes = column_filter(header_cols)
    else:
        keep_indexes = list(range(len(header_cols))) if header_cols else []

    if header_cols:
        filtered_header_cols = [header_cols[idx] for idx in keep_indexes]
        filtered_header = "\t".join(filtered_header_cols)
        rows_to_sort = [[r[idx] for idx in keep_indexes] for r in data_rows]
    else:
        filtered_header = ""
        rows_to_sort = data_rows

    # Determine the column index to sort rows by (e.g. gene_id or column index)
    if isinstance(sort_key_col, int):
        sort_idx = sort_key_col if sort_key_col < len(rows_to_sort[0]) else 0
    elif header_cols and sort_key_col in filtered_header_cols:
        sort_idx = filtered_header_cols.index(sort_key_col)
    else:
        sort_idx = 0 if rows_to_sort else None

    # Sort rows alphabetically to eliminate order-of-writing nondeterminism
    if sort_idx is not None:
        rows_sorted = sorted(rows_to_sort, key=lambda r: r[sort_idx])
    else:
        rows_sorted = rows_to_sort

    normalized_rows = ["\t".join(row) for row in rows_sorted]

    # Return standard unified TSV content string
    if filtered_header:
        return "\n".join([filtered_header] + normalized_rows).strip()
    return "\n".join(normalized_rows).strip()


# Normalizes Pyhmmer's space-delimited tblout files.
# It splits on whitespace at most 18 times (yielding exactly 19 fields) to preserve
# trailing annotation text columns that may contain spaces, and sorts rows by the first column.
def normalize_tblout(path: Path) -> str:
    return normalize_tsv(
        path,
        sort_key_col=0,
        line_parser=lambda line: line.split(maxsplit=18),
        has_header=False,
    )


# Normalizes the parsed annotation files by sorting rows alphabetically by 'gene_id'.
def normalize_parsed(path: Path) -> str:
    return normalize_tsv(path, sort_key_col="gene_id")


# Normalizes final annotation TSV files to make them comparable across run modes.
# Because certain databases are disabled in minimal or pyhmmer-only testing modes,
# their respective columns will be absent in actual outputs. This filter drops any database-specific
# columns that are NOT currently enabled, while keeping general headers and enabled database columns.
def normalize_annotation(
    path: Path,
    enabled_db_columns: set[str],
    all_db_columns: set[str],
) -> str:
    def col_filter(header_cols):
        return [
            idx
            for idx, col in enumerate(header_cols)
            if col in enabled_db_columns or col not in all_db_columns
        ]

    return normalize_tsv(path, sort_key_col="gene_id", column_filter=col_filter)


# Compares an actual generated output file to an expected fixture file.
# Supports optional normalizers (such as normalize_tsv or normalize_tblout) to prevent
# sorting or dynamic formatting issues from failing assertions.
# If files mismatch, it prints a unified diff showing the first 10 changes.
def compare_files(actual: Path, expected: Path, label: str, normalizer=None) -> bool:
    if not expected.exists():
        print(f"      {Colors.RED}Expected file not found: {expected.name}{Colors.END}")
        return False
    if not actual.exists():
        print(f"      {Colors.RED}Actual file not found: {actual}{Colors.END}")
        return False

    try:
        # Load and normalize content as required
        if normalizer is not None:
            actual_text = normalizer(actual)
            expected_text = normalizer(expected)
        else:
            actual_text = actual.read_text(encoding="utf-8", errors="replace")
            expected_text = expected.read_text(encoding="utf-8", errors="replace")
    except Exception as exc:
        print(f"      {Colors.RED}{exc}{Colors.END}")
        return False

    if actual_text == expected_text:
        return True

    # Generate and display a unified diff if the text mismatches
    diff = difflib.unified_diff(
        expected_text.splitlines(),
        actual_text.splitlines(),
        fromfile=str(expected),
        tofile=str(actual),
        lineterm="",
    )
    print(f"      {Colors.YELLOW}Content mismatch in {expected.name}:{Colors.END}")
    diff_lines = list(diff)
    # Show first 10 diff lines
    for line in diff_lines[:10]:
        print(f"        {line}")
    if len(diff_lines) > 10:
        print(f"        ... ({len(diff_lines) - 10} more differences)")
    return False


# Recursively compares all files in an actual output directory with an expected fixture directory.
# Identifies missing expected files, extra unexpected files, and mismatches within common files
# using a custom comparator or a standard line comparison.
def compare_directories(
    expected_root: Path,
    actual_root: Path,
    label: str,
    file_comparator=None,
) -> bool:
    if not expected_root.exists():
        print(
            f"      {Colors.RED}Expected directory not found: {expected_root}{Colors.END}"
        )
        return False
    if not actual_root.exists():
        print(
            f"      {Colors.RED}Actual output directory not found: {actual_root}{Colors.END}"
        )
        return False

    # Retrieve all files recursively, keeping paths relative to directory root
    expected_files = {
        path.relative_to(expected_root)
        for path in expected_root.rglob("*")
        if path.is_file()
    }
    actual_files = {
        path.relative_to(actual_root)
        for path in actual_root.rglob("*")
        if path.is_file()
    }
    success = True

    # Identify structural mismatches in directory contents
    missing = sorted(expected_files - actual_files)
    extra = sorted(actual_files - expected_files)
    if missing:
        print(f"      {Colors.RED}Missing expected files:{Colors.END}")
        for path in missing[:5]:
            print(f"        - {path}")
        if len(missing) > 5:
            print(f"        ... and {len(missing) - 5} more")
        success = False
    if extra:
        print(f"      {Colors.YELLOW}Unexpected extra files:{Colors.END}")
        for path in extra[:5]:
            print(f"        - {path}")
        if len(extra) > 5:
            print(f"        ... and {len(extra) - 5} more")
        success = False

    # Compare file contents for all files present in both directories
    for rel_path in sorted(expected_files & actual_files):
        expected_file = expected_root / rel_path
        actual_file = actual_root / rel_path
        if file_comparator is not None:
            if not file_comparator(actual_file, expected_file, f"{label}:{rel_path}"):
                success = False
        else:
            if not compare_files(actual_file, expected_file, f"{label}:{rel_path}"):
                success = False

    return success


# Standard wrapper around compare_directories for default directory verification.
def compare_directory(expected_root: Path, actual_root: Path, label: str) -> bool:
    return compare_directories(expected_root, actual_root, label)


# A specialized directory comparator for the final annotation directory.
# Uses a custom normalizer for TSV files to filter out disabled database columns
# dynamically, preventing false failures when comparing across different test modes.
def compare_annotation_directory(
    expected_root: Path,
    actual_root: Path,
    enabled_db_columns: set[str],
    all_db_columns: set[str],
    label: str,
) -> bool:
    def annotation_comparator(actual: Path, expected: Path, file_label: str) -> bool:
        # Use annotation normalizer for TSVs to ignore columns from disabled databases
        if expected.suffix == ".tsv":

            def normalizer(path):
                return normalize_annotation(path, enabled_db_columns, all_db_columns)

            return compare_files(actual, expected, file_label, normalizer=normalizer)
        return compare_files(actual, expected, file_label)

    return compare_directories(
        expected_root, actual_root, label, file_comparator=annotation_comparator
    )


# Cleans and normalizes nucleotide or amino acid sequences by removing newlines,
# whitespace, and converting all alphabetic characters to uppercase.
def normalize_sequence(raw: str) -> str:
    return "".join(ch for ch in raw if ch.isalpha()).upper()


# Extracts the amino acid /translation sequence from a GenBank (GBK) file
# for a given gene identifier (locus_tag).
def extract_gbk_translation(path: Path, locus_tag: str) -> str:
    if not path.exists():
        raise FileNotFoundError(f"Missing GBK: {path}")
    lines = path.read_text().splitlines()
    for idx, line in enumerate(lines):
        # Locate the specific locus tag entry
        if f'/locus_tag="{locus_tag}"' not in line:
            continue
        # Scan forward to locate the translation attribute
        for j in range(idx, len(lines)):
            if '/translation="' not in lines[j]:
                continue
            fragment = lines[j].split('/translation="', 1)[1]
            parts = [fragment]
            # If sequence is single-line, extract and return
            if '"' in fragment:
                return normalize_sequence(fragment.split('"', 1)[0])
            # For multi-line sequences, keep reading until the closing quote
            for k in range(j + 1, len(lines)):
                chunk = lines[k].strip()
                if '"' in chunk:
                    parts.append(chunk.split('"', 1)[0])
                    return normalize_sequence("".join(parts))
                parts.append(chunk)
            break
        break
    return ""


# Generates a short, human-readable preview of a protein sequence for logging.
def format_sequence_preview(sequence: str) -> str:
    preview = sequence[:20]
    suffix = "..." if len(sequence) > 20 else ""
    return f"{len(sequence)}aa:{preview}{suffix}"


# Compares protein sequences between a genome-based FAA and a protein-based FAA file.
# Rather than asserting identical line order, this matches sequence frequency counts,
# making comparisons order-invariant and showing clear mismatch statistics upon failure.
def compare_faa_sequences(genome_faa: Path, protein_faa: Path, label: str) -> bool:
    try:
        genome_seqs = list(parse_fasta(genome_faa).values())
        protein_seqs = list(parse_fasta(protein_faa).values())
    except FileNotFoundError as exc:
        print(f"      {Colors.RED}{exc}{Colors.END}")
        return False

    # Sequence ordering can be non-deterministic, so we count frequency of unique sequences
    genome_counts = {}
    for seq in genome_seqs:
        genome_counts[seq] = genome_counts.get(seq, 0) + 1

    protein_counts = {}
    for seq in protein_seqs:
        protein_counts[seq] = protein_counts.get(seq, 0) + 1

    if genome_counts == protein_counts:
        return True

    # Identify and display exact sequence differences
    missing = []
    for seq, count in genome_counts.items():
        diff = count - protein_counts.get(seq, 0)
        if diff > 0:
            missing.extend([seq] * diff)

    extra = []
    for seq, count in protein_counts.items():
        diff = count - genome_counts.get(seq, 0)
        if diff > 0:
            extra.extend([seq] * diff)

    print(f"      {Colors.YELLOW}Sequence mismatch:{Colors.END}")
    print(f"        Genome-only sequences: {len(missing)}")
    print(f"        Protein-only sequences: {len(extra)}")
    if missing:
        previews = ", ".join(format_sequence_preview(seq) for seq in missing[:3])
        print(f"        Genome-only previews: {previews}")
    if extra:
        previews = ", ".join(format_sequence_preview(seq) for seq in extra[:3])
        print(f"        Protein-only previews: {previews}")
    return False


# Parses a GFF3 file to extract tabular feature rows.
# Ignores headers, comments, and trailing inline FASTA blocks. Checks that every
# feature line contains exactly 9 columns as required by standard GFF3 specs.
def parse_gff_features(path: Path) -> list[list[str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing GFF: {path}")
    features = []
    for line in path.read_text().splitlines():
        # Stop parsing if we reach the FASTA section
        if line.startswith("##FASTA") or line.startswith(">"):
            break
        # Skip empty lines and comment lines
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            raise ValueError(f"Invalid GFF row with {len(fields)} columns: {path}")
        features.append(fields)
    return features


# Performs coordinate integrity validation on GFF feature entries.
# Ensures that indices are valid 1-based integers and that the start position
# is less than or equal to the end position.
def check_gff_coordinates(path: Path) -> bool:
    try:
        features = parse_gff_features(path)
    except (FileNotFoundError, ValueError) as exc:
        print(f"      {Colors.RED}{exc}{Colors.END}")
        return False

    for fields in features:
        try:
            start = int(fields[3])
            end = int(fields[4])
        except ValueError:
            print(
                f"      {Colors.RED}Non-integer coordinates in {path}: {fields[3]}-{fields[4]}{Colors.END}"
            )
            return False
        # GFF coordinates are 1-based, start must be <= end
        if start < 1 or end < 1 or end < start:
            print(
                f"      {Colors.RED}Invalid coordinates in {path}: {start}-{end}{Colors.END}"
            )
            return False

    return True


# Validates strand encoding values in GFF files.
# Standard GFF features should have strands encoded as '+', '-', '.', or '?'.
# Structural CRISPR repeat region features are expected to be unstranded ('.').
def check_gff_strand_encoding(path: Path) -> bool:
    try:
        features = parse_gff_features(path)
    except (FileNotFoundError, ValueError) as exc:
        print(f"      {Colors.RED}{exc}{Colors.END}")
        return False

    for fields in features:
        strand = fields[6]
        feature_type = fields[2]
        # Repeat regions (like CRISPR arrays) should have unstranded dot ('.') strand values
        if feature_type == "repeat_region":
            if strand != ".":
                print(
                    f"      {Colors.RED}CRISPR strand not '.' in {path}: {strand}{Colors.END}"
                )
                return False
        elif strand not in ("+", "-", ".", "?"):
            print(f"      {Colors.RED}Invalid strand in {path}: {strand}{Colors.END}")
            return False

    return True


# Helper function that loads a TSV file and parses its contents as a list
# of dictionaries mapped by header column names.
def parse_tsv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing TSV: {path}")
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != len(header):
            raise ValueError(f"Invalid TSV row in {path}")
        rows.append(dict(zip(header, fields, strict=True)))
    return rows


# Verifies CRISPR array predictions, ensuring valid repeat counts (>=2),
# length consistency (end - start + 1), and that the internal array geometry
# (repeats * repeat_length + spacers * spacer_length) matches the total length.
def check_crispr_arrays(path: Path) -> bool:
    try:
        rows = parse_tsv_rows(path)
    except (FileNotFoundError, ValueError) as exc:
        print(f"      {Colors.RED}{exc}{Colors.END}")
        return False

    if not rows:
        return True

    for row in rows:
        try:
            start = int(row["start"])
            end = int(row["end"])
            length = int(row["length"])
            repeat_length = int(row["repeat_length"])
            spacer_length = int(row["spacer_length"])
            num_repeats = int(row["num_repeats"])
        except (KeyError, ValueError):
            print(
                f"      {Colors.RED}Invalid CRISPR numeric fields in {path}{Colors.END}"
            )
            return False

        # Biologically valid CRISPR arrays require at least 2 repeats
        if num_repeats < 2:
            print(
                f"      {Colors.RED}CRISPR repeats < 2 in {path}: {row.get('crispr_id', '?')}{Colors.END}"
            )
            return False

        # Total predicted length should match the coordinate boundaries
        expected_length = end - start + 1
        if length != expected_length:
            print(
                f"      {Colors.RED}CRISPR length mismatch in {path}: {row.get('crispr_id', '?')}{Colors.END}"
            )
            return False

        # Geometry consistency: (num_repeats * repeat_length) + (num_spacers * spacer_length)
        inferred_length = (num_repeats * repeat_length) + (
            (num_repeats - 1) * spacer_length
        )
        if abs(expected_length - inferred_length) > (num_repeats - 1):
            print(
                f"      {Colors.RED}CRISPR geometry mismatch in {path}: {row.get('crispr_id', '?')}{Colors.END}"
            )
            return False

    return True


def check_trna_records(path: Path) -> bool:
    """
    Validate tRNA records found in the assembled `*_rna.tsv` file.

    Pyrococcus furiosus (DSM 3638) includes a singular loci for
    tRNA-His (Gene ID: 41713857) and tRNA-Cys (Gene ID: 41713313),
    respectively. We check for the exact location and sequences.
    """
    try:
        rows = parse_tsv_rows(path)
    except (FileNotFoundError, ValueError) as exc:
        print(f"      {Colors.RED}{exc}{Colors.END}")
        return False

    # Extract all rows corresponding to tRNA features
    trna_rows = [row for row in rows if row.get("rna_type", "").startswith("tRNA-")]
    if not trna_rows:
        return True

    success = True
    # Validate general structural properties of all parsed tRNAs
    for row in trna_rows:
        rna_type = row.get("rna_type", "")
        anti_codon = row.get("anti_codon", "")

        # tRNA identifiers should follow the standard pattern 'tRNA-[AminoAcidThreeLetters]'
        if not rna_type.startswith("tRNA-") or len(rna_type) < 6:
            print(
                f"      {Colors.RED}Invalid tRNA type in {path}: {rna_type}{Colors.END}"
            )
            return False
        # Verify proper upper-casing of the amino acid abbreviation
        if not rna_type[5].isupper():
            print(
                f"      {Colors.RED}Invalid tRNA type casing in {path}: {rna_type}{Colors.END}"
            )
            return False
        # Anticodon should either be undefined ('.') or a 3-letter lowercase code
        if anti_codon != ".":
            if anti_codon != anti_codon.lower() or len(anti_codon) != 3:
                print(
                    f"      {Colors.RED}Invalid tRNA anticodon in {path}: {anti_codon}{Colors.END}"
                )
                return False

    # Target reference tRNAs with known published coordinates/sequences for validation
    expected_trnas = [
        {
            "rna_type": "tRNA-His",
            "start": "1844790",
            "end": "1844866",
            "strand": "+",
            "sequence": "GCCGGGGTGGTGTAGCCTGGTTAGCACAGGGGACTGTGGATCCCCTGGCCCGGGTTCAAATCCCGGCCCCGGCCCCA",
            "reference": "Gene ID: 41713857",
        },
        {
            "rna_type": "tRNA-Cys",
            "start": "1380455",
            "end": "1380529",
            "strand": "-",
            "sequence": "GCCGGGATAGCCTAGAGGCCAGGCGGGGGACTGCAGATCCCCTTTACCCGGGTTCAAATCCCGGTCCCGGCTCCA",
            "reference": "Gene ID: 41713313",
        },
    ]

    # Assert that all expected reference tRNAs exist with exact attributes
    for expected in expected_trnas:
        matches = [
            row
            for row in trna_rows
            if row.get("rna_type") == expected["rna_type"]
            and row.get("start") == expected["start"]
            and row.get("end") == expected["end"]
            and row.get("strand") == expected["strand"]
            and row.get("sequence") == expected["sequence"]
        ]
        if not matches:
            print(
                f"      {Colors.RED}Missing or mismatched {expected['rna_type']} at "
                f"{expected['start']}-{expected['end']} ({expected['strand']}), "
                f"{expected['reference']} in {path}{Colors.END}"
            )
            success = False

    return success


def check_rrna_records(path: Path) -> bool:
    """
    Validate 5S rRNA records found in the assembled `*_rna.tsv` file.

    Pyrococcus furiosus (DSM 3638) encodes two 5S ribosomal RNA genes in its genome.
    We check for the exact location and sequence of the expected 5S rRNA loci.
    See https://www.ncbi.nlm.nih.gov/datasets/gene/41713297/ and
    https://www.ncbi.nlm.nih.gov/datasets/gene/41713504/ for details.
    """
    try:
        rows = parse_tsv_rows(path)
    except (FileNotFoundError, ValueError) as exc:
        print(f"      {Colors.RED}{exc}{Colors.END}")
        return False

    # Find 5S rRNA records - expect exactly two genes in Pyrococcus furiosus
    rrna_rows = [row for row in rows if row.get("rna_type", "") == "5S_rRNA"]
    if len(rrna_rows) != 2:
        print(
            f"      {Colors.RED}Expected 2 5S rRNA records in {path}, found {len(rrna_rows)}{Colors.END}"
        )
        return False

    # Define exact attributes of published 5S ribosomal rRNA coordinates and sequences
    expected_rrnas = [
        {
            "start": "1370591",
            "end": "1370712",
            "strand": "+",
            "sequence": "TACGGCGGCCATAGCGGGGGGGCCACACCCGGTCTCATTTCGAACCCGGAAGTTAAGCCCCCCAGCGATCCCGGTTGTACTGCCCTCCGAGAGGGGGCGGGAAGCCGGGGACGCCGCCGGCC",
        },
        {
            "start": "1541199",
            "end": "1541319",
            "strand": "-",
            "sequence": "TACGGCGGCCATAGCGGGGGGCCACACCCGGTCTCATTTCGAACCCGGAAGTTAAGCCCCCCAGCGATCCCGGTTGTACTGCCCTCCGAGAGGGGGCGGGAAGCCGGGGACGCCGCCGGCC",
        },
    ]

    success = True
    # Verify that each predicted 5S rRNA matches its expected genomic location and sequence exactly
    for expected in expected_rrnas:
        matches = [
            r
            for r in rrna_rows
            if r.get("start") == expected["start"]
            and r.get("end") == expected["end"]
            and r.get("strand") == expected["strand"]
            and r.get("sequence") == expected["sequence"]
        ]
        if not matches:
            print(
                f"      {Colors.RED}Missing or mismatched 5S rRNA at {expected['start']}-{expected['end']} ({expected['strand']}) in {path}{Colors.END}"
            )
            success = False

    return success


def check_argonaute_records(
    gbk_path: Path,
    gff_path: Path,
    tsv_path: Path,
    annotation_path: Path,
    fna_path: Path,
    faa_path: Path,
) -> bool:
    """
    Validate the published Argonaute gene record in Pyrococcus furiosus (DSM 3638).

    The PfAgo gene (Gene ID: 41712340) is reported on the minus strand at
    535964-538276 with a 2313 bp CDS. We check the feature coordinates, lengths,
    and sequences across GFF, GBK, TSV, and FASTA outputs.
    Reference: https://www.ncbi.nlm.nih.gov/datasets/gene/41712340/
    """
    feature_id = "JLPIBLAO_1_00524"
    expected_start = "535964"
    expected_end = "538276"
    expected_strand = "-"
    expected_length = "2313"

    # Exact expected nucleotide sequence of PfAgo (Gene ID: 41712340)
    expected_nucleotide = normalize_sequence(
        """
        ATGAAAGCGAAAGTTGTTATTAATCTGGTAAAGATAAATAAAAAAATAATTCCAGATAAAATATATGTAT
        ATAGGCTCTTTAATGACCCAGAGGAAGAACTTCAAAAAGAGGGATATTCTATTTATCGCCTTGCATATGA
        AAACGTTGGTATTGTTATTGATCCAGAGAATCTGATTATTGCCACCACAAAAGAGCTTGAATATGAGGGA
        GAATTTATTCCCGAGGGTGAAATATCATTCTCGGAATTGCGTAATGACTATCAGAGTAAACTTGTCTTAA
        GATTGTTGAAGGAGAATGGAATAGGTGAATATGAATTATCGAAACTTCTTAGGAAGTTTAGAAAACCAAA
        AACTTTTGGGGATTATAAGGTCATTCCAAGCGTTGAGATGAGTGTGATAAAGCATGATGAAGATTTTTAT
        TTAGTAATTCATATAATTCATCAAATACAATCAATGAAAACACTGTGGGAGCTTGTAAACAAAGACCCTA
        AAGAGCTTGAAGAGTTTTTAATGACTCACAAGGAAAATTTAATGCTGAAGGACATTGCGTCACCACTTAA
        AACTGTTTATAAGCCTTGCTTTGAAGAGTATACTAAAAAACCAAAACTTGATCATAATCAAGAAATAGTA
        AAGTACTGGTATAATTATCATATTGAAAGATATTGGAATACTCCAGAAGCTAAATTAGAGTTTTATAGAA
        AATTTGGCCAAGTTGATTTAAAACAACCTGCAATTCTAGCGAAATTTGCATCTAAAATAAAGAAAAACAA
        AAATTACAAGATATATCTCTTGCCCCAATTAGTTGTCCCAACTTATAACGCCGAACAATTAGAAAGTGAC
        GTGGCGAAAGAAATTTTAGAATATACAAAACTGATGCCAGAAGAACGTAAAGAGCTTTTAGAAAACATCC
        TTGCAGAAGTTGATAGTGACATTATAGACAAATCATTAAGTGAAATTGAAGTAGAAAAAATTGCACAAGA
        ATTGGAAAACAAAATAAGGGTTAGAGATGACAAAGGGAACAGTGTACCAATTTCTCAATTAAATGTTCAA
        AAATCTCAATTACTACTTTGGACGAATTATTCACGAAAATATCCAGTTATATTGCCTTATGAAGTTCCAG
        AAAAGTTCAGAAAAATACGGGAGATACCAATGTTTATAATTCTAGATTCAGGACTTCTTGCCGACATTCA
        AAATTTTGCAACTAATGAATTCAGGGAGTTAGTAAAGAGCATGTACTATAGTCTTGCCAAAAAATACAAT
        AGTCTTGCCAAAAAAGCAAGGTCAACAAATGAAATAGGATTACCCTTTTTAGACTTTCGTGGGAAAGAAA
        AGGTTATAACTGAAGATTTAAACTCTGACAAAGGTATTATAGAAGTTGTGGAGCAGGTATCTAGTTTCAT
        GAAAGGAAAAGAACTAGGCCTAGCCTTTATAGCTGCAAGAAATAAACTTTCATCTGAGAAGTTTGAGGAA
        ATTAAAAGGAGACTCTTTAATTTAAATGTAATATCTCAAGTGGTTAACGAAGATACTTTAAAAAATAAAA
        GAGACAAATACGATAGAAATCGACTTGATCTATTTGTCAGACACAACTTACTTTTTCAAGTATTATCTAA
        ACTTGGAGTAAAGTATTATGTGCTAGATTACAGGTTCAATTATGACTACATCATTGGAATTGATGTTGCT
        CCCATGAAGCGTTCTGAGGGATATATAGGTGGTAGTGCCGTAATGTTTGACTCTCAAGGTTACATACGAA
        AAATTGTCCCAATTAAAATTGGTGAACAAAGAGGAGAATCAGTTGACATGAATGAGTTTTTCAAAGAAAT
        GGTTGATAAATTTAAAGAGTTCAATATTAAGCTAGATAATAAGAAAATACTACTTCTAAGAGATGGAAGA
        ATTACTAATAACGAAGAAGAAGGGCTCAAGTATATTTCCGAGATGTTTGATATTGAGGTGGTCACGATGG
        ACGTGATAAAAAATCACCCTGTGAGAGCCTTTGCTAATATGAAAATGTACTTTAACCTAGGTGGTGCAAT
        ATACCTAATCCCCCATAAGCTCAAACAAGCTAAAGGAACACCAATACCAATAAAGTTAGCCAAAAAGAGG
        ATAATTAAGAATGGAAAAGTCGAGAAGCAGAGTATCACTAGACAGGACGTCTTAGATATCTTCATTTTAA
        CTCGTCTTAACTATGGGAGTATTTCTGCTGATATGCGGTTGCCAGCACCCGTTCACTATGCTCACAAGTT
        TGCAAATGCAATCAGAAATGAATGGAAAATAAAAGAAGAGTTCCTTGCTGAGGGATTTTTGTATTTTGTT
        TGA
        """
    )

    # Exact expected translation sequence of PfAgo protein (Gene ID: 41712340)
    expected_protein = normalize_sequence(
        """
        MKAKVVINLVKINKKIIPDKIYVYRLFNDPEEELQKEGYSIYRLAYENVGIVIDPENLIIATTKELEYEG
        EFIPEGEISFSELRNDYQSKLVLRLLKENGIGEYELSKLLRKFRKPKTFGDYKVIPSVEMSVIKHDEDFY
        LVIHIIHQIQSMKTLWELVNKDPKELEEFLMTHKENLMLKDIASPLKTVYKPCFEEYTKKPKLDHNQEIV
        KYWYNYHIERYWNTPEAKLEFYRKFGQVDLKQPAILAKFASKIKKNKNYKIYLLPQLVVPTYNAEQLESD
        VAKEILEYTKLMPEERKELLENILAEVDSDIIDKSLSEIEVEKIAQELENKIRVRDDKGNSVPISQLNVQ
        KSQLLLWTNYSRKYPVILPYEVPEKFRKIREIPMFIILDSGLLADIQNFATNEFRELVKSMYYSLAKKYN
        SLAKKARSTNEIGLPFLDFRGKEKVITEDLNSDKGIIEVVEQVSSFMKGKELGLAFIAARNKLSSEKFEE
        IKRRLFNLNVISQVVNEDTLKNKRDKYDRNRLDLFVRHNLLFQVLSKLGVKYYVLDYRFNYDYIIGIDVA
        PMKRSEGYIGGSAVMFDSQGYIRKIVPIKIGEQRGESVDMNEFFKEMVDKFKEFNIKLDNKKILLLRDGR
        ITNNEEEGLKYISEMFDIEVVTMDVIKNHPVRAFANMKMYFNLGGAIYLIPHKLKQAKGTPIPIKLAKKR
        IIKNGKVEKQSITRQDVLDIFILTRLNYGSISADMRLPAPVHYAHKFANAIRNEWKIKEEFLAEGFLYFV
        """
    )

    success = True

    # 1. Verify GFF3 CDS record exist with matching coordinates, strand, and ID
    try:
        gff_features = parse_gff_features(gff_path)
    except (FileNotFoundError, ValueError) as exc:
        print(f"      {Colors.RED}{exc}{Colors.END}")
        return False

    gff_matches = [
        row
        for row in gff_features
        if row[2] == "CDS"
        and row[3] == expected_start
        and row[4] == expected_end
        and row[6] == expected_strand
        and f"ID={feature_id}" in row[8]
    ]
    if not gff_matches:
        print(
            f"      {Colors.RED}Missing Argonaute CDS in {gff_path}: {feature_id} {expected_start}-{expected_end} ({expected_strand}){Colors.END}"
        )
        success = False

    # 2. Verify parsed gene entry matches coordinates and length in standard TSV features file
    try:
        tsv_rows = parse_tsv_rows(tsv_path)
    except (FileNotFoundError, ValueError) as exc:
        print(f"      {Colors.RED}{exc}{Colors.END}")
        return False

    tsv_matches = [
        row
        for row in tsv_rows
        if row.get("gene_id") == feature_id
        and row.get("start") == expected_start
        and row.get("end") == expected_end
        and row.get("strand") == expected_strand
        and row.get("gene_length") == expected_length
    ]
    if not tsv_matches:
        print(
            f"      {Colors.RED}Missing Argonaute gene in {tsv_path}: {feature_id} length {expected_length}{Colors.END}"
        )
        success = False

    # 3. Verify target database annotation hit in final annotation TSV (arCOG03890 for test_pyhmmer)
    try:
        annotation_rows = parse_tsv_rows(annotation_path)
    except (FileNotFoundError, ValueError) as exc:
        print(f"      {Colors.RED}{exc}{Colors.END}")
        return False

    annotation_matches = [
        row
        for row in annotation_rows
        if row.get("gene_id") == feature_id
        and row.get("gene_length") == expected_length
        and row.get("test_pyhmmer_hit") == "arCOG03890"
    ]
    if not annotation_matches:
        print(
            f"      {Colors.RED}Missing Argonaute annotation hit in {annotation_path}: {feature_id} arCOG03890{Colors.END}"
        )
        success = False

    # 4. Verify coordinates and locus tag match in GenBank (.gbk) file format
    if not gbk_path.exists():
        print(f"      {Colors.RED}Missing GBK: {gbk_path}{Colors.END}")
        return False

    gbk_text = gbk_path.read_text()
    if f"complement({expected_start}..{expected_end})" not in gbk_text:
        print(
            f"      {Colors.RED}Missing Argonaute GBK coordinates in {gbk_path}: complement({expected_start}..{expected_end}){Colors.END}"
        )
        success = False
    if f'/locus_tag="{feature_id}"' not in gbk_text:
        print(
            f"      {Colors.RED}Missing Argonaute locus_tag in {gbk_path}: {feature_id}{Colors.END}"
        )
        success = False

    # 5. Extract translation sequence from GBK and compare with expected protein
    translation = extract_gbk_translation(gbk_path, feature_id)
    if not translation or translation != expected_protein:
        print(
            f"      {Colors.RED}Argonaute GBK translation mismatch in {gbk_path}{Colors.END}"
        )
        success = False

    # 6. Verify nucleotide sequence identity in FNA FASTA output
    try:
        fna_records = parse_fasta(fna_path)
    except FileNotFoundError as exc:
        print(f"      {Colors.RED}{exc}{Colors.END}")
        return False

    if normalize_sequence(fna_records.get(feature_id, "")) != expected_nucleotide:
        print(
            f"      {Colors.RED}Argonaute FNA sequence mismatch in {fna_path}{Colors.END}"
        )
        success = False

    # 7. Verify protein sequence identity in FAA FASTA output
    try:
        faa_records = parse_fasta(faa_path)
    except FileNotFoundError as exc:
        print(f"      {Colors.RED}{exc}{Colors.END}")
        return False

    if normalize_sequence(faa_records.get(feature_id, "")) != expected_protein:
        print(
            f"      {Colors.RED}Argonaute FAA sequence mismatch in {faa_path}{Colors.END}"
        )
        success = False

    return success


# The main orchestrator that sets up and runs the entire Prokanota test suite.
# It parses inputs, runs biological invariants unit tests, sets up database selections,
# invokes Snakemake, and performs rigorous integration checks on all outputs.
def main() -> int:
    args = parse_args()

    if args.very_verbose:
        args.verbose = True

    # 1. Run biological truths and pipeline constraint unit tests
    if not run_pytest():
        return 1

    # In 'minimal' mode, we only run the Pytest suite and exit early
    if args.mode == "minimal":
        return 0

    # 2. Select databases based on chosen testing mode
    if args.mode == "full":
        databases_path = DATABASES_FULL
    else:
        databases_path = DATABASES_PYHMMER

    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)

    selected_databases = load_databases(databases_path)
    enabled_databases = get_enabled_databases(selected_databases)
    enabled_db_names = [
        db.get("name") for db in enabled_databases if isinstance(db.get("name"), str)
    ]
    enabled_db_columns_map = get_database_columns(enabled_databases)
    enabled_db_columns = set(
        col for cols in enabled_db_columns_map.values() for col in cols
    )

    full_databases = load_databases(DATABASES_FULL)
    all_db_columns = get_all_db_columns(full_databases)

    # 3. Create config and invoke the Snakemake workflow
    config_path = build_snakemake_config(CONFIGFILE, databases_path, args.mode)
    run_snakemake(config_path, args.verbose, args.very_verbose)

    success = True
    total_checks = 0
    passed_checks = 0

    # 4. Perform directory and file comparisons for each sample
    for sample_id in SAMPLES:
        print(f"\n{Colors.BOLD}{Colors.BLUE}Sample: {sample_id}{Colors.END}")

        # Check each database's raw and parsed output
        if not enabled_db_names:
            print(f"{Colors.YELLOW}  No databases enabled for this sample{Colors.END}")
            continue

        print(f"  Databases enabled: {', '.join(enabled_db_names)}")
        for db_name in enabled_db_names:
            print(f"\n  {Colors.BOLD}{db_name}:{Colors.END}")

            # Check raw hits output (tblout or tsv)
            hits_name = get_hits_filename(db_name)
            hits_normalizer = get_hits_normalizer(db_name)
            actual_tblout = INTERIM_DIR / sample_id / db_name / hits_name
            expected_tblout = EXPECTED_BASE_DIR / sample_id / db_name / hits_name
            total_checks += 1
            raw_ok = compare_files(
                actual_tblout,
                expected_tblout,
                f"raw:{sample_id}:{db_name}",
                normalizer=hits_normalizer,
            )
            if raw_ok:
                passed_checks += 1
                print(f"    {Colors.GREEN}✓ Raw hits{Colors.END}")
            else:
                success = False
                print(f"    {Colors.RED}✗ Raw hits{Colors.END}")

            # Check parsed database annotations
            actual_parsed = INTERIM_DIR / sample_id / db_name / "parsed_annotation.tsv"
            expected_parsed = (
                EXPECTED_BASE_DIR / sample_id / db_name / "parsed_annotation.tsv"
            )
            total_checks += 1
            if compare_files(
                actual_parsed,
                expected_parsed,
                f"parsed:{sample_id}:{db_name}",
                normalizer=normalize_parsed,
            ):
                passed_checks += 1
                print(f"    {Colors.GREEN}✓ Parsed annotation{Colors.END}")
            else:
                success = False
                print(f"    {Colors.RED}✗ Parsed annotation{Colors.END}")

        # Check final output feature files (GFF, GBK, TSV, FNA, FAA)
        print(f"\n  {Colors.BOLD}Final outputs:{Colors.END}")
        expected_features = EXPECTED_BASE_DIR / sample_id / "features"
        actual_features = RESULTS_DIR / sample_id / "features"
        total_checks += 1
        if compare_directory(
            expected_features, actual_features, f"features:{sample_id}"
        ):
            passed_checks += 1
            print(f"    {Colors.GREEN}✓ Features directory{Colors.END}")
        else:
            success = False
            print(f"    {Colors.RED}✗ Features directory{Colors.END}")

        # Check final annotated tables
        expected_annotation = EXPECTED_BASE_DIR / sample_id / "annotation"
        actual_annotation = RESULTS_DIR / sample_id / "annotation"
        total_checks += 1
        if compare_annotation_directory(
            expected_annotation,
            actual_annotation,
            enabled_db_columns,
            all_db_columns,
            f"annotation:{sample_id}",
        ):
            passed_checks += 1
            print(f"    {Colors.GREEN}✓ Annotation directory{Colors.END}")
        else:
            success = False
            print(f"    {Colors.RED}✗ Annotation directory{Colors.END}")

    # 5. Master table validation: Verify common table contains all sample rows
    print(f"\n{Colors.BOLD}{Colors.BLUE}Master table checks:{Colors.END}")
    common_annotation = RESULTS_DIR / "common" / "annotation" / "finalAnnotation.tsv"
    sample_annotations = [
        RESULTS_DIR / sample_id / "annotation" / f"{sample_id}_finalAnnotation.tsv"
        for sample_id in SAMPLES
    ]
    total_checks += 1
    try:
        common_rows = len(parse_tsv_rows(common_annotation))
        sample_rows = sum(len(parse_tsv_rows(path)) for path in sample_annotations)
        if common_rows == sample_rows:
            passed_checks += 1
            print(f"  {Colors.GREEN}✓ Common annotation row count{Colors.END}")
        else:
            success = False
            print(
                f"  {Colors.RED}✗ Common annotation row count "
                f"({common_rows} != {sample_rows}){Colors.END}"
            )
    except (FileNotFoundError, ValueError) as exc:
        success = False
        print(f"  {Colors.RED}✗ Common annotation row count: {exc}{Colors.END}")

    # 6. Cross-sample validation: Verify protein FASTA identity
    print(f"\n{Colors.BOLD}{Colors.BLUE}Cross-sample checks:{Colors.END}")
    genome_id = SAMPLES[0]
    protein_id = SAMPLES[1]
    genome_faa = RESULTS_DIR / genome_id / "features" / f"{genome_id}.faa"
    protein_faa = RESULTS_DIR / protein_id / "features" / f"{protein_id}.faa"
    total_checks += 1
    if compare_faa_sequences(genome_faa, protein_faa, "faa:genome-vs-protein"):
        passed_checks += 1
        print(f"  {Colors.GREEN}✓ Genome vs protein FAA sequences{Colors.END}")
    else:
        success = False
        print(f"  {Colors.RED}✗ Genome vs protein FAA sequences{Colors.END}")

    # 7. GFF structure validations (integer indices and correct strand encodings)
    print(f"\n{Colors.BOLD}{Colors.BLUE}Structural validation:{Colors.END}")
    for sample_id in SAMPLES:
        print(f"  {sample_id}:")
        sample_gff_path = RESULTS_DIR / sample_id / "features" / f"{sample_id}.gff"

        total_checks += 1
        if check_gff_coordinates(sample_gff_path):
            passed_checks += 1
            print(f"    {Colors.GREEN}✓ GFF coordinates{Colors.END}")
        else:
            success = False
            print(f"    {Colors.RED}✗ GFF coordinates{Colors.END}")

        total_checks += 1
        if check_gff_strand_encoding(sample_gff_path):
            passed_checks += 1
            print(f"    {Colors.GREEN}✓ GFF strand encoding{Colors.END}")
        else:
            success = False
            print(f"    {Colors.RED}✗ GFF strand encoding{Colors.END}")

    print(f"\n{Colors.BOLD}{Colors.BLUE}Biological validation:{Colors.END}")

    # 8. Verify specialized features and biological truths (CRISPR geometry, tRNAs, rRNAs, PfAgo)
    print(f"\n  {genome_sample} (specialized):")
    total_checks += 1
    if check_crispr_arrays(crispr_path):
        passed_checks += 1
        print(f"    {Colors.GREEN}✓ CRISPR arrays{Colors.END}")
    else:
        success = False
        print(f"    {Colors.RED}✗ CRISPR arrays{Colors.END}")

    total_checks += 1
    if check_trna_records(rna_path):
        passed_checks += 1
        print(f"    {Colors.GREEN}✓ tRNA records{Colors.END}")
    else:
        success = False
        print(f"    {Colors.RED}✗ tRNA records{Colors.END}")

    total_checks += 1
    if check_rrna_records(rna_path):
        passed_checks += 1
        print(f"    {Colors.GREEN}✓ rRNA records{Colors.END}")
    else:
        success = False
        print(f"    {Colors.RED}✗ rRNA records{Colors.END}")

    total_checks += 1
    if check_argonaute_records(
        gbk_path, gff_path, tsv_path, annotation_path, fna_path, faa_path
    ):
        passed_checks += 1
        print(f"    {Colors.GREEN}✓ Argonaute records{Colors.END}")
    else:
        success = False
        print(f"    {Colors.RED}✗ Argonaute records{Colors.END}")

    # 9. Report final pass/fail summary
    print(f"\n{Colors.BOLD}{Colors.BLUE}{'=' * 60}{Colors.END}")
    if success:
        print(
            f"{Colors.GREEN}{Colors.BOLD}✓ ALL TESTS PASSED ({passed_checks}/{total_checks} checks){Colors.END}"
        )
    else:
        print(
            f"{Colors.RED}{Colors.BOLD}✗ SOME TESTS FAILED ({passed_checks}/{total_checks} checks passed){Colors.END}"
        )
    print(f"{Colors.BOLD}{Colors.BLUE}{'=' * 60}{Colors.END}")

    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
