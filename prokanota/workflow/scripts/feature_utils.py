"""
Shared helper utilities for feature generation and protein metrics.
"""

import hashlib
from collections.abc import Iterable
from pathlib import Path


def generate_prokanota_id(
    sample_id: str,
    mode: str,
    sequences: Iterable[str],
) -> str:
    """
    Generate Prokanota's ten-letter ID for a genome or protein collection.

    The ID is based on the following information, in this exact order:

    1. The user-supplied sample ID.
    2. The input mode: ``"genome"`` or ``"protein"``.
    3. Every sequence, in the same order as in the input FASTA file.

    Before hashing, sequences are converted to uppercase. A single terminal
    ``*`` is removed from protein sequences because it marks the end of
    translation and is not an amino acid in the protein written by Prokanota.

    The length of each normalized sequence is added immediately before the
    sequence itself. This preserves sequence boundaries: two contigs ``AC``
    and ``GT`` must not be treated like the two contigs ``ACG`` and ``T``,
    even though both pairs would produce ``ACGT`` if simply concatenated.

    SHA-256 converts all of this information into a digest. The digest is then
    represented by exactly ten uppercase letters from A to Z. Consequently,
    the letters do not encode biological information; together they are a
    compact, deterministic representation of the hashed input.

    Args:
        sample_id: User-supplied sample identifier.
        mode: Either ``"genome"`` or ``"protein"``.
        sequences: Genome contigs or protein sequences in FASTA input order.

    Returns:
        A deterministic identifier containing exactly ten uppercase letters.

    Raises:
        ValueError: If mode is not ``"genome"`` or ``"protein"``.
    """
    if mode not in {"genome", "protein"}:
        raise ValueError("mode must be either 'genome' or 'protein'.")

    # SHA-256 processes a stream of bytes. Adding these values in separate
    # steps is equivalent to hashing sample_id + mode + all sequence fields.
    hash_object = hashlib.sha256()
    hash_object.update(sample_id.encode("utf-8"))
    hash_object.update(mode.encode("ascii"))

    for sequence in sequences:
        # FASTA sequences are case-insensitive, so upper- and lowercase input
        # must generate the same Prokanota ID.
        normalized_sequence = sequence.strip().upper()

        # A terminal '*' is a notation for the translation stop, not a protein
        # residue. Prokanota removes it from imported protein output, so it is
        # also removed before the ID is calculated. Internal '*' characters
        # are retained because they indicate a different input sequence.
        if mode == "protein" and normalized_sequence.endswith("*"):
            normalized_sequence = normalized_sequence[:-1]

        sequence_bytes = normalized_sequence.encode("ascii")

        # Store the sequence length as an eight-byte integer before storing
        # the sequence. Eight bytes allow any realistic sequence length. The
        # byte order is fixed so that all computers hash the same bytes.
        sequence_length = len(sequence_bytes).to_bytes(8, byteorder="big")
        hash_object.update(sequence_length)
        hash_object.update(sequence_bytes)

    # A SHA-256 digest is a large binary number. Taking its remainder after
    # division by 26^10 selects one value from all possible ten-letter IDs.
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    identifier_length = 10
    identifier_number = int.from_bytes(hash_object.digest(), byteorder="big") % (
        len(alphabet) ** identifier_length
    )

    # Convert that number to base 26. This is analogous to converting a number
    # to decimal, except that the "digits" are A through Z instead of 0 to 9.
    # divmod() returns both the remaining number and the next letter position.
    identifier_characters = []
    for _ in range(identifier_length):
        identifier_number, letter_index = divmod(identifier_number, len(alphabet))
        identifier_characters.append(alphabet[letter_index])

    # Base conversion finds the final letter first, so reverse the collected
    # letters to obtain the conventional left-to-right representation.
    identifier_characters.reverse()
    return "".join(identifier_characters)


def protein_molecular_weight(seq: str) -> float:
    """
    Estimate the molecular weight of a protein sequence in Daltons.
    Uses average residue masses (in Da) and adds the mass of one water molecule.
    This method is described in GASTEIGER, Elisabeth, et al. The proteomics protocols handbook, 2005, S. 571-607.
    The values can be found at https://web.expasy.org/findmod/findmod_masses.html#AA
    """
    weights = {
        "A": 71.0788,
        "R": 156.1875,
        "N": 114.1038,
        "D": 115.0886,
        "C": 103.1388,
        "E": 129.1155,
        "Q": 128.1307,
        "G": 57.0519,
        "H": 137.1411,
        "I": 113.1594,
        "L": 113.1594,
        "K": 128.1741,
        "M": 131.1926,
        "F": 147.1766,
        "P": 97.1167,
        "S": 87.0782,
        "T": 101.1051,
        "W": 186.2132,
        "Y": 163.1760,
        "V": 99.1326,
        "U": 150.0388,  # Selenocysteine
        "O": 237.3018,  # Pyrrolysine
    }
    if not seq:
        return 0.0
    mw = 0.0
    for aa in seq:
        if aa not in weights:
            return 0.0  # Unknown amino acid yields a zero weight result.
        mw += weights[aa]
    # Add the weight of one water molecule (approximately 18.015 Da) for the free termini
    mw += 18.01528
    return mw


def ensure_empty_file(path: str) -> None:
    """
    Create a truly empty file at the given path, ensuring parent directories exist.

    Args:
        path (str): File path to create.
    """
    file_path = Path(path)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    file_path.write_bytes(b"")


def reverse_complement(seq: str) -> str:
    """Reverse complement the given DNA sequence.

    This implementation is IUPAC-aware (see https://www.bioinformatics.org/sms/iupac.html)
    and preserves the case of the input sequence. Unknown characters are preserved as-is.
    """
    complement = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "N": "N",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "V": "B",
        "D": "H",
        "H": "D",
    }

    rc = []
    for ch in reversed(seq):
        up = ch.upper()
        comp = complement.get(up, up)
        rc.append(comp.lower() if ch.islower() else comp)
    return "".join(rc)


def get_sequence(contig_seq: str, start: int, end: int, strand: str) -> str:
    """Extract a sequence from a contig and reverse-complement if needed."""
    seq = contig_seq[start - 1 : end]
    if strand == "-":
        seq = reverse_complement(seq)
    return seq
