"""Biologically grounded tests for published truths."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
PACKAGE_ROOT = REPO_ROOT / "prokanota"
SCRIPTS_DIR = PACKAGE_ROOT / "workflow" / "scripts"

for path in (str(REPO_ROOT), str(SCRIPTS_DIR)):
    if path not in sys.path:
        sys.path.insert(0, path)

from prokanota.workflow.scripts import feature_utils

def _column_to_list(df, column):
    series = df[column]
    if hasattr(series, "to_list"):
        return series.to_list()
    return series.tolist()

# ---
# 1. Protein Molecular Weight Calculation (feature_utils.py)
# ---

# Test molecular weight calculation against ExPASy reference values.
# The protein_molecular_weight function uses average isotopic masses from:
# https://web.expasy.org/findmod/findmod_masses.html
# MW = sum(residue masses) + 18.01528 Da (one water molecule for free termini)
# Reference: Gasteiger et al., The Proteomics Protocols Handbook, 2005, pp. 571-607

def test_single_amino_acid_weights():
    """
    Verify individual amino acid residue weights match ExPASy reference.
    Each single amino acid's MW = residue mass + water (18.01528 Da).
    """
    # Glycine (smallest AA) -> residue mass 57.0519 Da (ExPASy)
    # Expected: 57.0519 + 18.01528 = 75.06718 Da
    assert abs(feature_utils.protein_molecular_weight("G") - 75.07) < 0.01

    # Tryptophan (largest AA) -> residue mass 186.2132 Da (ExPASy)
    # Expected: 186.2132 + 18.01528 = 204.22848 Da
    assert abs(feature_utils.protein_molecular_weight("W") - 204.23) < 0.01

def test_pfu_ago_protein_molecular_weight():
    """
    Validate MW calculation against Pyrococcus furiosus Argonaute (PfAgo).
    Sequence from https://www.uniprot.org/uniprotkb/Q8U3D2/entry
    Uniprot lists MW as 90390 Da, ExPASy ProtParam calculates 90390.15 Da.
    """
    pf_ago = "MKAKVVINLVKINKKIIPDKIYVYRLFNDPEEELQKEGYSIYRLAYENVGIVIDPENLIIATTKELEYEGEFIPEGEISFSELRNDYQSKLVLRLLKENGIGEYELSKLLRKFRKPKTFGDYKVIPSVEMSVIKHDEDFYLVIHIIHQIQSMKTLWELVNKDPKELEEFLMTHKENLMLKDIASPLKTVYKPCFEEYTKKPKLDHNQEIVKYWYNYHIERYWNTPEAKLEFYRKFGQVDLKQPAILAKFASKIKKNKNYKIYLLPQLVVPTYNAEQLESDVAKEILEYTKLMPEERKELLENILAEVDSDIIDKSLSEIEVEKIAQELENKIRVRDDKGNSVPISQLNVQKSQLLLWTNYSRKYPVILPYEVPEKFRKIREIPMFIILDSGLLADIQNFATNEFRELVKSMYYSLAKKYNSLAKKARSTNEIGLPFLDFRGKEKVITEDLNSDKGIIEVVEQVSSFMKGKELGLAFIAARNKLSSEKFEEIKRRLFNLNVISQVVNEDTLKNKRDKYDRNRLDLFVRHNLLFQVLSKLGVKYYVLDYRFNYDYIIGIDVAPMKRSEGYIGGSAVMFDSQGYIRKIVPIKIGEQRGESVDMNEFFKEMVDKFKEFNIKLDNKKILLLRDGRITNNEEEGLKYISEMFDIEVVTMDVIKNHPVRAFANMKMYFNLGGAIYLIPHKLKQAKGTPIPIKLAKKRIIKNGKVEKQSITRQDVLDIFILTRLNYGSISADMRLPAPVHYAHKFANAIRNEWKIKEEFLAEGFLYFV"
    mw = feature_utils.protein_molecular_weight(pf_ago)
    assert abs(mw - 90390) < 1.0 # allow for small rounding differences

def test_unknown_amino_acid_handling():
    """
    Verify that unknown amino acids (e.g., 'X', 'B', 'Z') return 0 Da.
    This alerts users to potential sequence issues rather than silently
    computing incorrect weights.
    """
    mw_with_unknown = feature_utils.protein_molecular_weight("MXK")
    assert mw_with_unknown == 0.0, "Unknown amino acid should yield zero weight"

def test_empty_sequence_molecular_weight():
    """
    Empty sequence should return zero.
    """
    assert feature_utils.protein_molecular_weight("") == 0

# ---
# 2. DNA Reverse Complement + Sequence Extraction (feature_utils.py)
# ---

# Test reverse complement function for biological correctness.
# Reverse complement is fundamental for handling minus-strand genes.
# The complement pairs are: A<->T, C<->G (Watson-Crick base pairing).

def test_reverse_complement_basic():
    """
    Verify basic reverse complement transformation.
    Input:  5'-ATGC-3'
    Output: 5'-GCAT-3'
    """
    assert feature_utils.reverse_complement("ATGC") == "GCAT"


def test_reverse_complement_palindrome():
    """
    Test palindromic sequence (restriction site).
    EcoRI recognition site (GAATTC) is palindromic - its reverse complement
    equals itself.
    """
    ecori_site = "GAATTC"
    assert feature_utils.reverse_complement(ecori_site) == ecori_site


def test_reverse_complement_preserves_length():
    """Ensure reverse complement preserves sequence length."""
    seq = "ATGCATGCATGC"
    assert len(feature_utils.reverse_complement(seq)) == len(seq)


def test_reverse_complement_uppercase_only():
    """
    Pyrodigal outputs uppercase sequences, so reverse_complement is designed to
    handle uppercase ACGT. Ambiguous bases (N, R, Y, etc.) are not handled at all and will just be reversed.
    Lowercase input will be reverse complemented but stay lowercase.
    """
    # Lowercase is reverse complemented but not converted to uppercase
    assert feature_utils.reverse_complement("atgc") == "gcat"  # NOT "GCAT"

    # Ambiguous bases pass through unchanged
    assert feature_utils.reverse_complement("ATNG") == "CNAT"  # N not complemented

def test_reverse_complement_empty():
    """Empty sequence should return empty string."""
    assert feature_utils.reverse_complement("") == ""

def test_single_nucleotide_reverse_complement():
    """Single nucleotide reverse complement is just the complement."""
    assert feature_utils.reverse_complement("A") == "T"
    assert feature_utils.reverse_complement("T") == "A"
    assert feature_utils.reverse_complement("G") == "C"
    assert feature_utils.reverse_complement("C") == "G"


def test_get_sequence_plus_strand():
    """Extract a sequence on the plus strand without reverse complement."""
    contig = "ATGCCGTA"
    assert feature_utils.get_sequence(contig, 2, 5, "+") == "TGCC"


def test_get_sequence_minus_strand():
    """Extract a sequence on the minus strand with reverse complement."""
    contig = "ATGCCGTA"
    assert feature_utils.get_sequence(contig, 2, 5, "-") == "GGCA"

# ---
# 3. E-value Filtering and Best Hit Selection (diamond_parse.py, pyhmmer_parse.py)
# ---

# Test e-value filtering logic for biological relevance.
# E-value represents the expected number of hits by chance; lower = more significant.
# Standard cutoffs: 1e-3 (permissive), 1e-5 (moderate), 1e-10 (stringent).
# Reference: Altschul et al., J Mol Biol 1990; 215(3):403-410 (BLAST paper)

def test_evalue_filtering_removes_insignificant_hits():
    """
    Verify that hits above e-value cutoff are filtered out.
    A hit with e-value 0.01 should be removed when cutoff is 1e-3.
    """
    pytest.importorskip("polars")
    from prokanota.workflow.scripts import diamond_parse
    import polars as pl

    hits_df = pl.DataFrame({
        "gene_id": ["gene1", "gene2", "gene3"],
        "accession": ["acc1", "acc2", "acc3"],
        "evalue": [1e-10, 1e-3, 0.01],  # significant, borderline, insignificant
        "score": [100.0, 50.0, 30.0],
    })

    filtered = diamond_parse.filter_and_deduplicate(hits_df, evalue_cutoff=1e-3)

    # Only gene1 and gene2 should remain (evalue <= 1e-3)
    assert len(filtered) == 2
    assert "gene3" not in filtered["gene_id"].to_list()


def test_best_hit_selection_by_score():
    """
    When a gene has multiple hits, select the one with highest bitscore.
    Bitscore is normalized and database-size independent, making it
    preferable for comparing hits across different databases.
    Reference: Karlin & Altschul, PNAS 1990; 87(6):2264-2268
    """
    pytest.importorskip("polars")
    from prokanota.workflow.scripts import diamond_parse
    import polars as pl

    hits_df = pl.DataFrame({
        "gene_id": ["gene1", "gene1", "gene1"],
        "accession": ["pfam1", "pfam2", "pfam3"],
        "evalue": [1e-20, 1e-25, 1e-15],
        "score": [150.0, 200.0, 100.0],  # pfam2 has highest score
    })

    filtered = diamond_parse.filter_and_deduplicate(hits_df, evalue_cutoff=1e-3)

    assert len(filtered) == 1
    assert filtered["accession"][0] == "pfam2"  # Highest score wins


def test_tiebreaker_uses_evalue():
    """
    When scores are tied, prefer the hit with lower (better) e-value.
    This ensures we select the most statistically significant hit.
    """
    pytest.importorskip("polars")
    from prokanota.workflow.scripts import diamond_parse
    import polars as pl

    hits_df = pl.DataFrame({
        "gene_id": ["gene1", "gene1"],
        "accession": ["hit_a", "hit_b"],
        "evalue": [1e-30, 1e-20],  # hit_a has better e-value
        "score": [100.0, 100.0],   # Same score
    })

    filtered = diamond_parse.filter_and_deduplicate(hits_df, evalue_cutoff=1e-3)

    assert filtered["accession"][0] == "hit_a"  # Lower e-value wins tie


# ---
# 4. Mapping File Validation (mapping_utils.py)
# ---

# Test mapping file parsing for robustness with real-world database formats.
# Mapping files link database accessions to functional annotations.
# Common issues: inconsistent empty markers, duplicate accessions, malformed rows.


def test_empty_marker_normalization():
    """
    Verify that various empty markers are normalized to canonical '*'.
    Different databases use different conventions for missing data:
    - Pfam: '-' for no clan
    - COG: 'NA' for unassigned
    - KEGG: empty string for no pathway
    All should be normalized to '*' for consistent downstream processing.
    """
    pytest.importorskip("pandas")
    from prokanota.workflow.scripts import mapping_utils
    content = (
        "ACC001\t-\tDescription\tCategory\n"
        "ACC002\tNA\tDescription\tCategory\n"
        "ACC003\t\tDescription\tCategory\n"
        "ACC004\tn/a\tDescription\tCategory\n"
    )

    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as handle:
            handle.write(content)
            tmp_path = Path(handle.name)

        df = mapping_utils.parse_mapping_file(tmp_path)

        # All empty markers should be normalized to '*'
        assert _column_to_list(df, "short_name") == ["*", "*", "*", "*"]
    finally:
        if tmp_path is not None:
            tmp_path.unlink(missing_ok=True)


def test_duplicate_accession_rejection():
    """
    Verify that duplicate accessions raise an error.
    Duplicate accessions would cause ambiguous annotation assignments
    and must be caught during validation, not silently ignored.
    """
    pytest.importorskip("pandas")
    from prokanota.workflow.scripts import mapping_utils
    content = (
        "PF00001\tGPCR\tG-protein coupled receptor\tSignaling\n"
        "PF00001\tGPCR_dup\tDuplicate entry\tSignaling\n"
    )

    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as handle:
            handle.write(content)
            tmp_path = Path(handle.name)

        with pytest.raises(ValueError, match="Duplicate accession"):
            mapping_utils.parse_mapping_file(tmp_path)
    finally:
        if tmp_path is not None:
            tmp_path.unlink(missing_ok=True)


def test_pfam_accession_format():
    """
    Test parsing of Pfam-style accessions (PFxxxxx.version).
    Reference: https://pfam.xfam.org/help#tabview=tab4
    """
    pytest.importorskip("pandas")
    from prokanota.workflow.scripts import mapping_utils
    content = (
        "PF00001.23\t7tm_1\tGPCR family\tSignaling\n"
        "PF00002.28\t7tm_2\tSecretin family\tSignaling\n"
    )

    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as handle:
            handle.write(content)
            tmp_path = Path(handle.name)

        df = mapping_utils.parse_mapping_file(tmp_path)

        assert _column_to_list(df, "accession") == ["PF00001.23", "PF00002.28"]
    finally:
        if tmp_path is not None:
            tmp_path.unlink(missing_ok=True)


def test_wrong_column_count_rejection():
    """
    Verify malformed rows with wrong column counts raise an error.
    """
    pytest.importorskip("pandas")
    from prokanota.workflow.scripts import mapping_utils
    content = (
        "ACC001\tGeneA\tDescA\tcatA\n"
        "ACC002\tGeneB\tDescB\n"
    )

    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as handle:
            handle.write(content)
            tmp_path = Path(handle.name)

        with pytest.raises(ValueError, match="Expected 4 columns, found 3"):
            mapping_utils.parse_mapping_file(tmp_path)
    finally:
        if tmp_path is not None:
            tmp_path.unlink(missing_ok=True)


def test_missing_accession_rejection():
    """
    Verify empty accession keys are rejected.
    """
    pytest.importorskip("pandas")
    from prokanota.workflow.scripts import mapping_utils
    content = (
        "\tGeneA\tDescA\tcatA\n"
        "ACC002\tGeneB\tDescB\tcatB\n"
    )

    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as handle:
            handle.write(content)
            tmp_path = Path(handle.name)

        with pytest.raises(ValueError, match=r"Column 1 \(accession\) must not be empty"):
            mapping_utils.parse_mapping_file(tmp_path)
    finally:
        if tmp_path is not None:
            tmp_path.unlink(missing_ok=True)


def test_wrapped_double_quotes_are_stripped():
    """
    Verify wrapped double quotes are stripped from all fields.
    """
    pytest.importorskip("pandas")
    from prokanota.workflow.scripts import mapping_utils
    content = '"CDD:XXXXX"\t"gene_1"\t"description with spaces"\t"category"\n'

    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as handle:
            handle.write(content)
            tmp_path = Path(handle.name)

        df = mapping_utils.parse_mapping_file(tmp_path)

        assert _column_to_list(df, "accession") == ["CDD:XXXXX"]
        assert _column_to_list(df, "short_name") == ["gene_1"]
        assert _column_to_list(df, "description") == ["description with spaces"]
        assert _column_to_list(df, "category") == ["category"]
    finally:
        if tmp_path is not None:
            tmp_path.unlink(missing_ok=True)


def test_special_characters_are_preserved():
    """
    Verify description contents with special characters are preserved.
    """
    pytest.importorskip("pandas")
    from prokanota.workflow.scripts import mapping_utils
    content = "ACC001\tGeneA\tcontains 'single' and \"double\" quotes Ω\tcatA\n"

    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as handle:
            handle.write(content)
            tmp_path = Path(handle.name)

        df = mapping_utils.parse_mapping_file(tmp_path)

        assert _column_to_list(df, "description") == ["contains 'single' and \"double\" quotes Ω"]
    finally:
        if tmp_path is not None:
            tmp_path.unlink(missing_ok=True)


def test_empty_markers_normalized_across_all_fields():
    """
    Verify empty markers are normalized across non-key fields.
    """
    pytest.importorskip("pandas")
    from prokanota.workflow.scripts import mapping_utils
    content = (
        "ACC001\t\tNA\tNULL\n"
        "ACC002\t-\tN/A\t*\n"
        "ACC003\t  \tnull\t   \n"
    )

    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as handle:
            handle.write(content)
            tmp_path = Path(handle.name)

        df = mapping_utils.parse_mapping_file(tmp_path)

        assert _column_to_list(df, "short_name") == ["*", "*", "*"]
        assert _column_to_list(df, "description") == ["*", "*", "*"]
        assert _column_to_list(df, "category") == ["*", "*", "*"]
    finally:
        if tmp_path is not None:
            tmp_path.unlink(missing_ok=True)


def test_validate_mapping_file_returns_stats():
    """
    Verify validation returns row count and sample accessions.
    """
    pytest.importorskip("pandas")
    from prokanota.workflow.scripts import mapping_utils
    content = (
        "ACC001\tGeneA\tDescA\tcatA\n"
        "ACC002\tGeneB\tDescB\tcatB\n"
        "ACC003\tGeneC\tDescC\tcatC\n"
    )

    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as handle:
            handle.write(content)
            tmp_path = Path(handle.name)

        stats = mapping_utils.validate_mapping_file(tmp_path)

        assert stats["row_count"] == 3
        assert stats["sample_accessions"] == ["ACC001", "ACC002", "ACC003"]
    finally:
        if tmp_path is not None:
            tmp_path.unlink(missing_ok=True)


# ---
# 5. Annotation Merging Priority (merge_annotations.py)
# ---

# Test annotation merging respects database priority order.
# When multiple databases annotate the same gene, priority determines
# which annotation appears first in the output.
# Higher priority (lower order number) databases should be processed first.

def test_merge_preserves_database_order():
    """
    Verify annotations are merged in specified priority order.
    Example: Pfam (order=1) should appear before COG (order=2).
    This ensures consistent, reproducible annotation output.
    """
    pytest.importorskip("polars")
    from prokanota.workflow.scripts import merge_annotations
    import polars as pl

    base_df = pl.DataFrame({
        "gene_id": ["gene1", "gene2"],
        "contig_id": ["c1", "c1"],
    })

    pfam_content = "gene_id\tpfam_hit\n" "gene1\tPF00001\n"
    cog_content = "gene_id\tcog_hit\n" "gene1\tCOG0001\n" "gene2\tCOG0002\n"

    pfam_path = None
    cog_path = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as handle:
            handle.write(pfam_content)
            pfam_path = Path(handle.name)
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as handle:
            handle.write(cog_content)
            cog_path = Path(handle.name)

        merged = merge_annotations.merge_annotations(
            base_df,
            [("pfam", pfam_path), ("cog", cog_path)],
            gene_id_column="gene_id",
        )

        columns = merged.columns
        # Verify both columns exist and are in correct order
        assert "pfam_hit" in columns, "pfam_hit column missing after merge"
        assert "cog_hit" in columns, "cog_hit column missing after merge"
        assert columns.index("pfam_hit") < columns.index("cog_hit")
    finally:
        if pfam_path is not None:
            pfam_path.unlink(missing_ok=True)
        if cog_path is not None:
            cog_path.unlink(missing_ok=True)


def test_merge_fills_missing_with_asterisk():
    """
    Verify genes without hits in a database get '*' placeholder.
    This maintains consistent column structure across all genes
    and clearly indicates absence of annotation vs. missing data.
    """
    pytest.importorskip("polars")
    from prokanota.workflow.scripts import merge_annotations
    import polars as pl

    base_df = pl.DataFrame({
        "gene_id": ["gene1", "gene2", "gene3"],
    })

    annotation_content = "gene_id\tdb_hit\n" "gene1\tHIT001\n"

    annotation_path = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as handle:
            handle.write(annotation_content)
            annotation_path = Path(handle.name)

        merged = merge_annotations.merge_annotations(
            base_df,
            [("db", annotation_path)],
            gene_id_column="gene_id",
        )

        assert merged.filter(pl.col("gene_id") == "gene2")["db_hit"][0] == "*"
        assert merged.filter(pl.col("gene_id") == "gene3")["db_hit"][0] == "*"
    finally:
        if annotation_path is not None:
            annotation_path.unlink(missing_ok=True)


# ---
# 6. ID Generation Determinism (feature_utils.py)
# ---

# Test that genome/gene IDs are deterministic and reproducible.
# MD5-based ID generation ensures same input always produces same output.
# This is critical for reproducibility and tracking samples across runs.


def test_genome_id_deterministic():
    """
    Verify same sample_id + sequence produces identical genome_id.
    Reproducibility is essential for comparing results across pipeline runs.
    """
    sample_id = "test_sample"

    hash1 = feature_utils.hash_sample_id(sample_id)
    hash2 = feature_utils.hash_sample_id(sample_id)

    assert hash1 == hash2

    id1 = feature_utils.map_hexdigest_to_id(hash1, length=8)
    id2 = feature_utils.map_hexdigest_to_id(hash2, length=8)

    assert id1 == id2


def test_genome_id_uniqueness():
    """
    Verify different sample_ids produce different genome_ids.
    MD5 collision resistance ensures unique identifiers.
    """
    id1 = feature_utils.map_hexdigest_to_id(
        feature_utils.hash_sample_id("sample_A"), length=8
    )
    id2 = feature_utils.map_hexdigest_to_id(
        feature_utils.hash_sample_id("sample_B"), length=8
    )

    assert id1 != id2


def test_hexdigest_mapping_alphabet():
    """
    Verify hex-to-alphabetic mapping produces only uppercase letters.
    The mapping converts 0-9 to A-J and a-f to K-P, ensuring IDs
    are filesystem-safe and visually distinct from accession numbers.
    """
    # MD5 produces lowercase hex: 0-9, a-f
    test_hex = "0123456789abcdef"
    mapped = feature_utils.map_hexdigest_to_id(test_hex, length=16)

    # Should produce: ABCDEFGHIJKLMNOP
    assert mapped == "ABCDEFGHIJKLMNOP"
    assert mapped.isalpha()
    assert mapped.isupper()
