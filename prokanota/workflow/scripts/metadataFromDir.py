#!/usr/bin/env python3

import argparse
import csv
import re
import sys
from pathlib import Path

FASTA_EXTENSIONS = [
    ".fasta",
    ".fa",
    ".fna",
    ".faa",
    ".ffn",
]


def normalize_filename(filename):
    """
    Normalize filename into a sample_id.

    Removes the FASTA extension and replaces spaces/special characters
    with underscores.
    """
    name = filename

    # Remove known FASTA extension
    for ext in FASTA_EXTENSIONS:
        if name.lower().endswith(ext):
            name = name[: -len(ext)]
            break

    # Replace any non-alphanumeric character with underscore
    normalized = re.sub(r"[^A-Za-z0-9_]+", "_", name)

    # Collapse repeated underscores
    normalized = re.sub(r"_+", "_", normalized)

    # Remove leading/trailing underscores
    normalized = normalized.strip("_")

    return normalized


def make_unique_sample_id(sample_id, used_sample_ids):
    """
    Ensure sample_id is unique.

    If sample_id already exists, append _2, _3, etc.
    """
    if sample_id not in used_sample_ids:
        used_sample_ids.add(sample_id)
        return sample_id

    counter = 2

    while True:
        candidate = f"{sample_id}_{counter}"

        if candidate not in used_sample_ids:
            used_sample_ids.add(candidate)
            return candidate

        counter += 1


def generate_metadata(input_dir, output_file, mode, comment="", force=False):
    """
    Generate metadata CSV from a directory of FASTA files.

    Args:
        input_dir: Path to directory containing FASTA files
        output_file: Path to output CSV file
        mode: Type of sequences ("dna" or "protein")
        comment: Optional comment for all entries
        force: Whether to overwrite existing output file

    Raises:
        ValueError: If input validation fails
    """
    input_dir = Path(input_dir)
    output_file = Path(output_file)

    # Check if input directory exists
    if not input_dir.is_dir():
        raise ValueError(f"Directory does not exist: {input_dir}")

    # Prevent accidental overwrite
    if output_file.exists() and not force:
        raise ValueError(
            f"Output file already exists: {output_file}. Use force=True to overwrite."
        )

    # Find all FASTA files
    fasta_files = []

    for path in input_dir.iterdir():
        if path.is_file() and any(
            path.name.lower().endswith(ext) for ext in FASTA_EXTENSIONS
        ):
            fasta_files.append(path)

    if not fasta_files:
        raise ValueError(f"No FASTA files found in: {input_dir}")

    fasta_files = sorted(fasta_files, key=lambda p: p.name.lower())

    rows = []
    used_sample_ids = set()

    for fasta_file in fasta_files:
        sample_id_base = normalize_filename(fasta_file.name)

        if not sample_id_base:
            raise ValueError(
                f"Could not generate a valid sample_id from: {fasta_file.name}"
            )

        sample_id = make_unique_sample_id(sample_id_base, used_sample_ids)

        # Always use absolute paths
        full_path = fasta_file.resolve()

        rows.append(
            [
                sample_id,
                str(full_path),
                mode,
                comment,
            ]
        )

    # Generate output file
    try:
        output_file.parent.mkdir(parents=True, exist_ok=True)

        with output_file.open("w", newline="") as outfile:
            writer = csv.writer(outfile)

            # Write header
            writer.writerow(["sample_id", "path", "input_type", "comment"])

            # Write data rows
            writer.writerows(rows)

        return f"Successfully created metadata file '{output_file}' with {len(rows)} entries."

    except OSError as e:
        raise ValueError(f"Could not write to output file '{output_file}': {e}") from e


def main():
    parser = argparse.ArgumentParser(
        description="Generate metadata CSV from directory of FASTA files"
    )

    parser.add_argument(
        "--path",
        required=True,
        type=Path,
        help="Directory containing FASTA files",
    )

    parser.add_argument(
        "--out",
        required=False,
        type=Path,
        default=None,
        help="Output CSV file path (default: metadata.csv in current directory)",
    )

    parser.add_argument(
        "--mode",
        required=True,
        choices=["dna", "protein"],
        help="Type of sequences: dna or protein",
    )

    parser.add_argument(
        "--comment",
        default="",
        help="Optional comment for all entries; default is empty string",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite output file if it already exists",
    )

    args = parser.parse_args()

    input_dir = args.path
    output_file = args.out if args.out is not None else Path.cwd() / "metadata.csv"

    try:
        message = generate_metadata(
            input_dir, output_file, args.mode, args.comment, args.force
        )
        print(message)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
