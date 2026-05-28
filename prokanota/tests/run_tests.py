#!/usr/bin/env python3
"""Run a minimal Snakemake test run and compare outputs to fixtures."""

from __future__ import annotations

import argparse
import difflib
import subprocess
import sys
from pathlib import Path
from collections import Counter

import yaml

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = PACKAGE_ROOT / "workflow" / "Snakefile"
CONFIGFILE = PACKAGE_ROOT / "tests" / "test-config.yaml"
DATABASES_FULL = PACKAGE_ROOT / "tests" / "test-databases.yaml"
DATABASES_PYHMMER = PACKAGE_ROOT / "tests" / "test-pyhmmer.yaml"
INTERIM_DIR = PACKAGE_ROOT / "tests" / "output" / "interim"
RESULTS_DIR = PACKAGE_ROOT / "tests" / "output" / "results"
EXPECTED_BASE_DIR = PACKAGE_ROOT / "tests" / "expected"
SAMPLES = [
    "Pyrococcus_furiosus_DSM_3638-GCF_008245085",
    "Pyrococcus_furiosus_DSM_3638-GCF_008245085_protein",
]

# ANSI color codes for terminal output
class Colors:
    GREEN = "\033[92m"
    RED = "\033[91m"
    YELLOW = "\033[93m"
    BLUE = "\033[94m"
    BOLD = "\033[1m"
    END = "\033[0m"


def get_hits_filename(db_name: str) -> str:
    if db_name == "test_pyhmmer":
        return "hits.tblout"
    return "hits.tsv"


def get_hits_normalizer(db_name: str):
    if db_name == "test_pyhmmer":
        return normalize_tblout
    return None


PYTEST_TESTS = PACKAGE_ROOT / "tests" / "test_published_biological_truths.py"


def load_yaml_file(path: Path, label: str) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Missing {label} file: {path}")
    data = yaml.safe_load(path.read_text())
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ValueError(f"Invalid {label} file {path}: root must be a mapping")
    return data


def load_databases(path: Path) -> list[dict[str, object]]:
    data = load_yaml_file(path, "databases")
    databases = data.get("databases")
    if not isinstance(databases, list):
        raise ValueError(f"Invalid databases file {path}: missing databases list")
    return databases


def get_enabled_databases(databases: list[dict[str, object]]) -> list[dict[str, object]]:
    enabled = []
    for db in databases:
        if not isinstance(db, dict):
            continue
        if not db.get("enabled", True):
            continue
        enabled.append(db)
    return enabled


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


def get_all_db_columns(databases: list[dict[str, object]]) -> set[str]:
    all_columns = set()
    for db in databases:
        for col in db.get("columns", []) if isinstance(db, dict) else []:
            if isinstance(col, dict):
                col_name = col.get("name")
                if isinstance(col_name, str):
                    all_columns.add(col_name)
    return all_columns


def build_snakemake_config(config_path: Path, databases_path: Path, mode_label: str) -> Path:
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


def run_pytest() -> bool:
    cmd = [sys.executable, "-m", "pytest", "-q", str(PYTEST_TESTS)]
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, cwd=PACKAGE_ROOT)
    return result.returncode == 0


def run_snakemake(configfile: Path, verbose: bool, very_verbose: bool) -> None:
    if very_verbose:
        verbose = True
    cmd = [
        "snakemake",
        "--snakefile",
        str(SNAKEFILE),
        "--configfile",
        str(configfile),
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


def normalize_annotation(
    path: Path,
    enabled_db_columns: set[str],
    all_db_columns: set[str],
) -> str:
    if not path.exists():
        raise FileNotFoundError(f"Missing annotation TSV: {path}")
    lines = [line.rstrip("\n") for line in path.read_text().splitlines() if line.strip()]
    if not lines:
        return ""
    header = lines[0]
    header_cols = header.split("\t")
    keep_indexes = [
        idx
        for idx, col in enumerate(header_cols)
        if col in enabled_db_columns or col not in all_db_columns
    ]
    filtered_header_cols = [header_cols[idx] for idx in keep_indexes]
    filtered_header = "\t".join(filtered_header_cols)

    rows = []
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != len(header_cols):
            raise ValueError(f"Invalid TSV row in {path}")
        rows.append([fields[idx] for idx in keep_indexes])

    if not rows:
        return filtered_header.strip()

    try:
        gene_idx = filtered_header_cols.index("gene_id")
    except ValueError:
        gene_idx = None

    if gene_idx is not None:
        rows_sorted = sorted(rows, key=lambda r: r[gene_idx])
    else:
        rows_sorted = rows

    normalized_rows = ["\t".join(row) for row in rows_sorted]
    return "\n".join([filtered_header] + normalized_rows).strip()


def compare_files(actual: Path, expected: Path, normalizer, label: str) -> bool:
    if not expected.exists():
        print(f"      {Colors.RED}Expected file not found: {expected.name}{Colors.END}")
        return False
    try:
        actual_text = normalizer(actual)
        expected_text = normalizer(expected)
    except (FileNotFoundError, ValueError) as exc:
        print(f"      {Colors.RED}{exc}{Colors.END}")
        return False

    if actual_text == expected_text:
        return True

    diff = difflib.unified_diff(
        expected_text.splitlines(),
        actual_text.splitlines(),
        fromfile=str(expected),
        tofile=str(actual),
        lineterm="",
    )
    print(f"      {Colors.YELLOW}Content mismatch (content differs in order or values):{Colors.END}")
    diff_lines = list(diff)
    # Show first 10 diff lines
    for line in diff_lines[:10]:
        print(f"        {line}")
    if len(diff_lines) > 10:
        print(f"        ... ({len(diff_lines) - 10} more differences)")
    return False


def compare_files_bytes(actual: Path, expected: Path, label: str) -> bool:
    if not expected.exists():
        print(f"      {Colors.RED}Expected file not found: {expected}{Colors.END}")
        return False
    if not actual.exists():
        print(f"      {Colors.RED}Actual file not found: {actual}{Colors.END}")
        return False

    actual_bytes = actual.read_bytes()
    expected_bytes = expected.read_bytes()
    if actual_bytes == expected_bytes:
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
    print(f"      {Colors.YELLOW}Byte mismatch in {expected.name}:{Colors.END}")
    diff_lines = list(diff)
    for line in diff_lines[:5]:
        print(f"        {line}")
    if len(diff_lines) > 5:
        print(f"        ... ({len(diff_lines) - 5} more lines differ)")
    return False


def compare_directory(expected_root: Path, actual_root: Path, label: str) -> bool:
    if not expected_root.exists():
        print(f"      {Colors.RED}Expected directory not found: {expected_root}{Colors.END}")
        return False
    if not actual_root.exists():
        print(f"      {Colors.RED}Actual output directory not found: {actual_root}{Colors.END}")
        return False

    expected_files = {path.relative_to(expected_root) for path in expected_root.rglob("*") if path.is_file()}
    actual_files = {path.relative_to(actual_root) for path in actual_root.rglob("*") if path.is_file()}
    success = True

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

    for rel_path in sorted(expected_files & actual_files):
        expected_file = expected_root / rel_path
        actual_file = actual_root / rel_path
        if not compare_files_bytes(actual_file, expected_file, f"{label}:{rel_path}"):
            success = False

    return success


def compare_annotation_directory(
    expected_root: Path,
    actual_root: Path,
    enabled_db_columns: set[str],
    all_db_columns: set[str],
    label: str,
) -> bool:
    if not expected_root.exists():
        print(f"      {Colors.RED}Expected directory not found: {expected_root}{Colors.END}")
        return False
    if not actual_root.exists():
        print(f"      {Colors.RED}Actual output directory not found: {actual_root}{Colors.END}")
        return False

    expected_files = {path.relative_to(expected_root) for path in expected_root.rglob("*") if path.is_file()}
    actual_files = {path.relative_to(actual_root) for path in actual_root.rglob("*") if path.is_file()}
    success = True

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

    for rel_path in sorted(expected_files & actual_files):
        expected_file = expected_root / rel_path
        actual_file = actual_root / rel_path
        if expected_file.suffix == ".tsv":
            normalizer = lambda path: normalize_annotation(path, enabled_db_columns, all_db_columns)
            if not compare_files(actual_file, expected_file, normalizer, f"{label}:{rel_path}"):
                success = False
        else:
            if not compare_files_bytes(actual_file, expected_file, f"{label}:{rel_path}"):
                success = False

    return success


def parse_faa_sequences(path: Path) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"Missing FASTA: {path}")
    sequences = []
    current = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if current:
                sequences.append("".join(current))
                current = []
            continue
        if line.strip():
            current.append(line.strip())
    if current:
        sequences.append("".join(current))
    return sequences


def format_sequence_preview(sequence: str) -> str:
    preview = sequence[:20]
    suffix = "..." if len(sequence) > 20 else ""
    return f"{len(sequence)}aa:{preview}{suffix}"


def compare_faa_sequences(genome_faa: Path, protein_faa: Path, label: str) -> bool:
    try:
        genome_seqs = parse_faa_sequences(genome_faa)
        protein_seqs = parse_faa_sequences(protein_faa)
    except FileNotFoundError as exc:
        print(f"      {Colors.RED}{exc}{Colors.END}")
        return False

    genome_counts = Counter(genome_seqs)
    protein_counts = Counter(protein_seqs)
    if genome_counts == protein_counts:
        return True

    missing = list((genome_counts - protein_counts).elements())
    extra = list((protein_counts - genome_counts).elements())
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


def parse_gff_features(path: Path) -> list[list[str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing GFF: {path}")
    features = []
    for line in path.read_text().splitlines():
        if line.startswith("##FASTA") or line.startswith(">"):
            break
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            raise ValueError(f"Invalid GFF row with {len(fields)} columns: {path}")
        features.append(fields)
    return features


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
            print(f"      {Colors.RED}Non-integer coordinates in {path}: {fields[3]}-{fields[4]}{Colors.END}")
            return False
        if start < 1 or end < 1 or end < start:
            print(f"      {Colors.RED}Invalid coordinates in {path}: {start}-{end}{Colors.END}")
            return False

    return True


def check_gff_strand_encoding(path: Path) -> bool:
    try:
        features = parse_gff_features(path)
    except (FileNotFoundError, ValueError) as exc:
        print(f"      {Colors.RED}{exc}{Colors.END}")
        return False

    for fields in features:
        strand = fields[6]
        feature_type = fields[2]
        if feature_type == "repeat_region":
            if strand != ".":
                print(f"      {Colors.RED}CRISPR strand not '.' in {path}: {strand}{Colors.END}")
                return False
        elif strand not in ("+", "-", ".", "?"):
            print(f"      {Colors.RED}Invalid strand in {path}: {strand}{Colors.END}")
            return False

    return True


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
        rows.append(dict(zip(header, fields)))
    return rows


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
            print(f"      {Colors.RED}Invalid CRISPR numeric fields in {path}{Colors.END}")
            return False

        if num_repeats < 2:
            print(f"      {Colors.RED}CRISPR repeats < 2 in {path}: {row.get('crispr_id', '?')}{Colors.END}")
            return False

        expected_length = (end - start + 1)
        if length != expected_length:
            print(f"      {Colors.RED}CRISPR length mismatch in {path}: {row.get('crispr_id', '?')}{Colors.END}")
            return False

        inferred_length = (num_repeats * repeat_length) + ((num_repeats - 1) * spacer_length)
        if abs(expected_length - inferred_length) > (num_repeats - 1):
            print(f"      {Colors.RED}CRISPR geometry mismatch in {path}: {row.get('crispr_id', '?')}{Colors.END}")
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

    trna_rows = [row for row in rows if row.get("rna_type", "").startswith("tRNA-")]
    if not trna_rows:
        return True

    success = True
    for row in trna_rows:
        rna_type = row.get("rna_type", "")
        anti_codon = row.get("anti_codon", "")

        if not rna_type.startswith("tRNA-") or len(rna_type) < 6:
            print(f"      {Colors.RED}Invalid tRNA type in {path}: {rna_type}{Colors.END}")
            return False
        if not rna_type[5].isupper():
            print(f"      {Colors.RED}Invalid tRNA type casing in {path}: {rna_type}{Colors.END}")
            return False
        if anti_codon != ".":
            if anti_codon != anti_codon.lower() or len(anti_codon) != 3:
                print(f"      {Colors.RED}Invalid tRNA anticodon in {path}: {anti_codon}{Colors.END}")
                return False

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

    # Find 5S rRNA records
    rrna_rows = [row for row in rows if row.get("rna_type", "") == "5S_rRNA"]
    if len(rrna_rows) != 2:
        print(f"      {Colors.RED}Expected 2 5S rRNA records in {path}, found {len(rrna_rows)}{Colors.END}")
        return False

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


def main() -> int:
    args = parse_args()

    if args.very_verbose:
        args.verbose = True

    if not run_pytest():
        return 1

    if args.mode == "minimal":
        return 0

    if args.mode == "full":
        databases_path = DATABASES_FULL
    else:
        databases_path = DATABASES_PYHMMER

    selected_databases = load_databases(databases_path)
    enabled_databases = get_enabled_databases(selected_databases)
    enabled_db_names = [db.get("name") for db in enabled_databases if isinstance(db.get("name"), str)]
    enabled_db_columns_map = get_database_columns(enabled_databases)
    enabled_db_columns = set(col for cols in enabled_db_columns_map.values() for col in cols)

    full_databases = load_databases(DATABASES_FULL)
    all_db_columns = get_all_db_columns(full_databases)

    config_path = build_snakemake_config(CONFIGFILE, databases_path, args.mode)
    run_snakemake(config_path, args.verbose, args.very_verbose)

    success = True
    total_checks = 0
    passed_checks = 0

    for sample_id in SAMPLES:
        print(f"\n{Colors.BOLD}{Colors.BLUE}Sample: {sample_id}{Colors.END}")
        
        # Check each database's raw and parsed output
        if not enabled_db_names:
            print(f"{Colors.YELLOW}  No databases enabled for this sample{Colors.END}")
            continue

        print(f"  Databases enabled: {', '.join(enabled_db_names)}")
        for db_name in enabled_db_names:
            print(f"\n  {Colors.BOLD}{db_name}:{Colors.END}")
            
            # Check raw hits output
            hits_name = get_hits_filename(db_name)
            hits_normalizer = get_hits_normalizer(db_name)
            actual_tblout = INTERIM_DIR / sample_id / db_name / hits_name
            expected_tblout = EXPECTED_BASE_DIR / sample_id / db_name / hits_name
            total_checks += 1
            if hits_normalizer:
                raw_ok = compare_files(actual_tblout, expected_tblout, hits_normalizer, f"raw:{sample_id}:{db_name}")
            else:
                raw_ok = compare_files_bytes(actual_tblout, expected_tblout, f"raw:{sample_id}:{db_name}")
            if raw_ok:
                passed_checks += 1
                print(f"    {Colors.GREEN}✓ Raw hits{Colors.END}")
            else:
                success = False
                print(f"    {Colors.RED}✗ Raw hits{Colors.END}")

            # Check parsed annotation
            actual_parsed = INTERIM_DIR / sample_id / db_name / "parsed_annotation.tsv"
            expected_parsed = EXPECTED_BASE_DIR / sample_id / db_name / "parsed_annotation.tsv"
            total_checks += 1
            if compare_files(actual_parsed, expected_parsed, normalize_parsed, f"parsed:{sample_id}:{db_name}"):
                passed_checks += 1
                print(f"    {Colors.GREEN}✓ Parsed annotation{Colors.END}")
            else:
                success = False
                print(f"    {Colors.RED}✗ Parsed annotation{Colors.END}")

        # Check features and annotation (these are the same for all databases)
        print(f"\n  {Colors.BOLD}Final outputs:{Colors.END}")
        expected_features = EXPECTED_BASE_DIR / sample_id / "features"
        actual_features = RESULTS_DIR / sample_id / "features"
        total_checks += 1
        if compare_directory(expected_features, actual_features, f"features:{sample_id}"):
            passed_checks += 1
            print(f"    {Colors.GREEN}✓ Features directory{Colors.END}")
        else:
            success = False
            print(f"    {Colors.RED}✗ Features directory{Colors.END}")

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

    # Genome vs protein comparison
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

    # GFF and structured data validation
    print(f"\n{Colors.BOLD}{Colors.BLUE}Structural validation:{Colors.END}")
    for sample_id in SAMPLES:
        print(f"  {sample_id}:")
        gff_path = RESULTS_DIR / sample_id / "features" / f"{sample_id}.gff"
        
        total_checks += 1
        if check_gff_coordinates(gff_path):
            passed_checks += 1
            print(f"    {Colors.GREEN}✓ GFF coordinates{Colors.END}")
        else:
            success = False
            print(f"    {Colors.RED}✗ GFF coordinates{Colors.END}")
        
        total_checks += 1
        if check_gff_strand_encoding(gff_path):
            passed_checks += 1
            print(f"    {Colors.GREEN}✓ GFF strand encoding{Colors.END}")
        else:
            success = False
            print(f"    {Colors.RED}✗ GFF strand encoding{Colors.END}")

    # Specialized feature checks (genome only)
    genome_sample = SAMPLES[0]
    crispr_path = RESULTS_DIR / genome_sample / "features" / f"{genome_sample}_crispr.tsv"
    rna_path = RESULTS_DIR / genome_sample / "features" / f"{genome_sample}_rna.tsv"
    
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

    # Summary
    print(f"\n{Colors.BOLD}{Colors.BLUE}{'='*60}{Colors.END}")
    if success:
        print(f"{Colors.GREEN}{Colors.BOLD}✓ ALL TESTS PASSED ({passed_checks}/{total_checks} checks){Colors.END}")
    else:
        print(f"{Colors.RED}{Colors.BOLD}✗ SOME TESTS FAILED ({passed_checks}/{total_checks} checks passed){Colors.END}")
    print(f"{Colors.BOLD}{Colors.BLUE}{'='*60}{Colors.END}")

    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
