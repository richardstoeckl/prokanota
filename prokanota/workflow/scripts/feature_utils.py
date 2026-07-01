"""
Shared helper utilities for feature generation and protein metrics.
"""

import hashlib
from pathlib import Path


def map_hexdigest_to_id(hexdigest: str, length: int = 8) -> str:
    """
    Map a hex digest to an alphabetic ID using the same mapping as features.py.

    Args:
        hexdigest (str): MD5 hex digest.
        length (int): Number of characters to use from the digest.

    Returns:
        str: Mapped ID string.
    """
    # Map 0-9 to A-J and a-f to K-P
    mapping = {str(i): chr(65 + i) for i in range(10)}
    mapping.update({chr(97 + i): chr(75 + i) for i in range(6)})
    return "".join(mapping[c] for c in hexdigest[:length])


def hash_sample_id(sample_id: str) -> str:
    """
    Hash the sample ID with MD5 and return the hex digest.

    Args:
        sample_id (str): The sample identifier.

    Returns:
        str: MD5 hex digest.
    """
    hash_obj = hashlib.md5(sample_id.encode("utf-8"))
    return hash_obj.hexdigest()


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
