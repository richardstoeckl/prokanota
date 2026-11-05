"""Shared text, FASTA and TSV I/O for workflow scripts."""

import csv
import gzip

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def open_text(path, mode="rt"):
    """Open UTF-8 text, using gzip compression for paths ending in .gz."""
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    # Preserve newline characters; the FASTA reader handles them explicitly,
    # and Python's CSV reader/writer needs newline="" for correct quoted fields.
    return opener(path, mode, encoding="utf-8", newline="")


def iter_fasta(path, *, allow_empty=False):
    """Read DNA or protein FASTA once, yielding one validated record at a time.

    Ignore blank lines and sequence-line edge whitespace, but reject internal
    whitespace. Join wrapped sequence lines and uppercase their contents.
    Checks cover file structure and ASCII sequences, not a DNA or protein alphabet;
    ambiguous residues, gaps and protein stops are not removed.

    Exhaust the iterator to validate the whole file. It keeps the current
    record and all seen IDs in memory, rather than every sequence. A caller
    that stops early must close the iterator to close the file.
    """
    with open_text(path) as handle:
        # Duplicate IDs make it ambiguous which sequence an annotation belongs to.
        seen = set()
        identifier = None
        title = ""
        header_line = 0
        fragments = []
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            # Never silently discard sequence text before the first header.
            # Check the original line so an indented '>' is not treated as valid.
            if identifier is None and not line.startswith(">"):
                raise ValueError(
                    f"FASTA {path}: line {line_number}: expected a header starting with '>'"
                )
            if line.startswith(">"):
                # A new header ends the previous record. Build that record now
                # instead of buffering the whole file or parsing it a second time.
                if identifier is not None:
                    yield _fasta_record(path, identifier, title, header_line, fragments)
                title = line[1:].rstrip()
                # Use the first word as the FASTA ID. Keep the title separately,
                # except trailing whitespace, so descriptions are preserved.
                words = title.split(maxsplit=1)
                identifier = words[0] if words else ""
                if not identifier:
                    raise ValueError(
                        f"FASTA {path}: line {line_number}: missing record ID"
                    )
                if identifier in seen:
                    raise ValueError(
                        f"FASTA {path}: line {line_number}: duplicate record ID {identifier!r}"
                    )
                seen.add(identifier)
                header_line = line_number
                fragments = []
            else:
                # An embedded '>' can become a false header when FASTA is rewritten.
                if ">" in line:
                    raise ValueError(
                        f"FASTA {path}: record {identifier!r}, line {line_number}: unexpected '>' in sequence"
                    )
                # Tolerate line wrapping and edge padding, but do not silently
                # repair internal spaces/tabs: that could hide a damaged sequence.
                sequence_line = line.strip()
                if any(char.isspace() for char in sequence_line):
                    raise ValueError(
                        f"FASTA {path}: record {identifier!r}, line {line_number}: whitespace inside sequence"
                    )
                # Check before uppercasing: some Unicode letters become ASCII
                # when uppercased, which would otherwise conceal invalid input.
                if not sequence_line.isascii():
                    raise ValueError(
                        f"FASTA {path}: record {identifier!r}, line {line_number}: non-ASCII sequence"
                    )
                fragments.append(sequence_line)
        if identifier is None:
            # Internal prediction output can legitimately contain no records.
            # This exception never permits a header with an empty sequence.
            if not allow_empty:
                raise ValueError(f"FASTA {path}: line 1: no records found")
            return
        # EOF, rather than another header, ends the last record.
        yield _fasta_record(path, identifier, title, header_line, fragments)


def _fasta_record(path, identifier, title, line_number, fragments):
    """Assemble one record and retain its header location for later diagnostics."""
    # Line wrapping and letter case should not change the sequence used downstream.
    # Keep '*' here: protein output and ID generation handle terminal stops
    # separately. Removing one here could make those callers remove a second one.
    sequence = "".join(fragments).upper()
    if not sequence:
        raise ValueError(
            f"FASTA {path}: record {identifier!r}, line {line_number}: empty sequence"
        )
    record = SeqRecord(Seq(sequence), id=identifier, name=identifier, description=title)
    # Later checks can report the source line without reopening the input file.
    record.annotations["fasta_line"] = line_number
    return record


def read_fasta(path, *, allow_empty=False):
    """Validate the entire FASTA and return its records as a reusable list."""
    # These callers need all records, or must detect malformed input before
    # writing output. Counters use iter_fasta directly to avoid storing sequences.
    return list(iter_fasta(path, allow_empty=allow_empty))


def check_source_proteins(path, records):
    """Reject source proteins that would be empty after terminal-stop removal."""
    # A single '*' passes structural FASTA checks but contains no amino acids.
    # Inspect the normalized value without modifying the records used for hashing.
    for record in records:
        if str(record.seq).removesuffix("*") == "":
            raise ValueError(
                f"FASTA {path}: record {record.id!r}, line {record.annotations['fasta_line']}: "
                "empty protein after removing terminal stop"
            )


def write_fasta_record(handle, header, sequence, width=None):
    """Write one FASTA record, optionally wrapping sequence lines at width."""
    # Callers choose the header and biological normalization. Serialization must
    # not change case or remove stops, which could alter the intended sequence.
    handle.write(f">{header}\n")
    sequence = str(sequence)
    if width:
        for i in range(0, len(sequence), width):
            handle.write(sequence[i : i + width] + "\n")
    else:
        handle.write(sequence + "\n")


def tsv_writer(handle):
    """Return a tab-separated writer that safely quotes tabs and newlines."""
    # FASTA descriptions can contain tabs or quotes. CSV quoting preserves those
    # characters inside one field instead of shifting the annotation columns.
    return csv.writer(handle, delimiter="\t", lineterminator="\n")
