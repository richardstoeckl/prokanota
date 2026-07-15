"""Biologically grounded tests for published truths."""

from __future__ import annotations

import importlib
import sys
import tempfile
from pathlib import Path

import pytest

from prokanota.workflow.scripts import feature_utils

REPO_ROOT = Path(__file__).resolve().parents[2]
PACKAGE_ROOT = REPO_ROOT / "prokanota"
SCRIPTS_DIR = PACKAGE_ROOT / "workflow" / "scripts"

for path in (str(REPO_ROOT), str(SCRIPTS_DIR)):
    if path not in sys.path:
        sys.path.insert(0, path)


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
    assert abs(mw - 90390) < 1.0  # allow for small rounding differences


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
    handle uppercase ACGT. Lowercase input will be reverse complemented but stay lowercase.
    """
    # Lowercase is reverse complemented but not converted to uppercase
    assert feature_utils.reverse_complement("atgc") == "gcat"  # NOT "GCAT"

    # Ambiguous bases
    assert feature_utils.reverse_complement("ATRC") == "GYAT"


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
    import polars as pl

    from prokanota.workflow.scripts import diamond_parse

    hits_df = pl.DataFrame(
        {
            "gene_id": ["gene1", "gene2", "gene3"],
            "accession": ["acc1", "acc2", "acc3"],
            "evalue": [1e-10, 1e-3, 0.01],  # significant, borderline, insignificant
            "score": [100.0, 50.0, 30.0],
            "sstart": [1, 1, 1],
            "send": [100, 100, 100],
            "subject_length": [100, 100, 100],
        }
    )

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
    import polars as pl

    from prokanota.workflow.scripts import diamond_parse

    hits_df = pl.DataFrame(
        {
            "gene_id": ["gene1", "gene1", "gene1"],
            "accession": ["pfam1", "pfam2", "pfam3"],
            "evalue": [1e-20, 1e-25, 1e-15],
            "score": [150.0, 200.0, 100.0],  # pfam2 has highest score
            "sstart": [1, 1, 1],
            "send": [100, 100, 100],
            "subject_length": [100, 100, 100],
        }
    )

    filtered = diamond_parse.filter_and_deduplicate(hits_df, evalue_cutoff=1e-3)

    assert len(filtered) == 1
    assert filtered["accession"][0] == "pfam2"  # Highest score wins


def test_tiebreaker_uses_evalue():
    """
    When scores are tied, prefer the hit with lower (better) e-value.
    This ensures we select the most statistically significant hit.
    """
    pytest.importorskip("polars")
    import polars as pl

    from prokanota.workflow.scripts import diamond_parse

    hits_df = pl.DataFrame(
        {
            "gene_id": ["gene1", "gene1"],
            "accession": ["hit_a", "hit_b"],
            "evalue": [1e-30, 1e-20],  # hit_a has better e-value
            "score": [100.0, 100.0],  # Same score
            "sstart": [1, 1],
            "send": [100, 100],
            "subject_length": [100, 100],
        }
    )

    filtered = diamond_parse.filter_and_deduplicate(hits_df, evalue_cutoff=1e-3)

    assert filtered["accession"][0] == "hit_a"  # Lower e-value wins tie


@pytest.mark.parametrize(
    "module_name",
    ["diamond_parse", "mmseqs2_parse", "rpsblast_parse", "pyhmmer_parse"],
)
def test_exact_tiebreakers_use_reference_coverage_then_identifier(module_name):
    """Exact score/E-value ties prefer coverage, then a stable identifier."""
    pytest.importorskip("polars")
    import polars as pl

    parser = importlib.import_module(f"prokanota.workflow.scripts.{module_name}")
    data = {
        "gene_id": ["coverage", "coverage", "identifier", "identifier"],
        "accession": ["hit_a", "hit_z", "hit_z", "hit_a"],
        "query_name": ["hit_a", "hit_z", "hit_z", "hit_a"],
        "evalue": [1e-20] * 4,
        "score": [100.0] * 4,
    }
    if module_name == "pyhmmer_parse":
        data["accession"] = ["-"] * 4
        data["reference_coverage"] = [0.5, 0.9, 0.9, 0.9]
    else:
        data.update(
            {
                "sstart": [1] * 4,
                "send": [50, 90, 90, 90],
                "subject_length": [100] * 4,
            }
        )

    filtered = parser.filter_and_deduplicate(pl.DataFrame(data))
    selected = dict(zip(filtered["gene_id"], filtered["query_name"], strict=True))

    assert selected == {"coverage": "hit_z", "identifier": "hit_a"}


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
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as handle:
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
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as handle:
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
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as handle:
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

    content = "ACC001\tGeneA\tDescA\tcatA\nACC002\tGeneB\tDescB\n"

    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as handle:
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

    content = "\tGeneA\tDescA\tcatA\nACC002\tGeneB\tDescB\tcatB\n"

    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as handle:
            handle.write(content)
            tmp_path = Path(handle.name)

        with pytest.raises(
            ValueError, match=r"Column 1 \(accession\) must not be empty"
        ):
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
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as handle:
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
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as handle:
            handle.write(content)
            tmp_path = Path(handle.name)

        df = mapping_utils.parse_mapping_file(tmp_path)

        assert _column_to_list(df, "description") == [
            "contains 'single' and \"double\" quotes Ω"
        ]
    finally:
        if tmp_path is not None:
            tmp_path.unlink(missing_ok=True)


def test_empty_markers_normalized_across_all_fields():
    """
    Verify empty markers are normalized across non-key fields.
    """
    pytest.importorskip("pandas")
    from prokanota.workflow.scripts import mapping_utils

    content = "ACC001\t\tNA\tNULL\nACC002\t-\tN/A\t*\nACC003\t  \tnull\t   \n"

    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as handle:
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
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as handle:
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
    import polars as pl

    from prokanota.workflow.scripts import merge_annotations

    base_df = pl.DataFrame(
        {
            "gene_id": ["gene1", "gene2"],
            "contig_id": ["c1", "c1"],
        }
    )

    pfam_content = "gene_id\tpfam_hit\ngene1\tPF00001\n"
    cog_content = "gene_id\tcog_hit\ngene1\tCOG0001\ngene2\tCOG0002\n"

    pfam_path = None
    cog_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as handle:
            handle.write(pfam_content)
            pfam_path = Path(handle.name)
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as handle:
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
    import polars as pl

    from prokanota.workflow.scripts import merge_annotations

    base_df = pl.DataFrame(
        {
            "gene_id": ["gene1", "gene2", "gene3"],
        }
    )

    annotation_content = "gene_id\tdb_hit\ngene1\tHIT001\n"

    annotation_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as handle:
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
# 6. Prokanota ID Generation (feature_utils.py)
# ---

# Test that genome and protein IDs are deterministic and reflect all biological
# input used to generate feature identifiers. This is critical for comparing
# annotations across runs without silently reusing an ID for changed input.


def test_prokanota_id_is_deterministic():
    """
    Verify the same sample, mode, and sequences always produce the same ID.
    Reproducibility is essential for comparing results across pipeline runs.
    """
    id1 = feature_utils.generate_prokanota_id(
        sample_id="test_sample",
        mode="genome",
        sequences=["ACGT", "TGCA"],
    )
    id2 = feature_utils.generate_prokanota_id(
        sample_id="test_sample",
        mode="genome",
        sequences=["ACGT", "TGCA"],
    )

    assert id1 == id2


def test_prokanota_id_known_result():
    """
    Verify a fixed reference example remains stable after implementation changes.
    This expected result fixes the exact SHA-256 input order, length encoding,
    and ten-letter conversion as part of Prokanota's reproducibility contract.
    """
    prokanota_id = feature_utils.generate_prokanota_id(
        sample_id="test_sample",
        mode="genome",
        sequences=["ACGT", "TGCA"],
    )

    assert prokanota_id == "IQPATDKZLQ"


def test_prokanota_id_uses_sample_id():
    """
    Verify changing the sample ID changes the generated Prokanota ID.
    The sample ID is intentionally part of the identifier calculation.
    """
    id1 = feature_utils.generate_prokanota_id(
        sample_id="sample_A",
        mode="genome",
        sequences=["ACGT"],
    )
    id2 = feature_utils.generate_prokanota_id(
        sample_id="sample_B",
        mode="genome",
        sequences=["ACGT"],
    )

    assert id1 != id2


def test_prokanota_id_contains_ten_uppercase_letters():
    """
    Verify Prokanota IDs contain exactly ten letters from A to Z.
    The letters are a compact representation of the SHA-256 hash and do not
    encode biological meaning individually.
    """
    prokanota_id = feature_utils.generate_prokanota_id(
        sample_id="test_sample",
        mode="genome",
        sequences=["ACGT"],
    )

    assert len(prokanota_id) == 10
    assert set(prokanota_id) <= set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")


def test_prokanota_id_preserves_sequence_boundaries():
    """
    Verify different contig boundaries produce different IDs.
    AC + GT and ACG + T both concatenate to ACGT, but they represent
    different genome assemblies and must remain distinguishable.
    """
    two_equal_contigs = feature_utils.generate_prokanota_id(
        sample_id="test_sample",
        mode="genome",
        sequences=["AC", "GT"],
    )
    unequal_contigs = feature_utils.generate_prokanota_id(
        sample_id="test_sample",
        mode="genome",
        sequences=["ACG", "T"],
    )

    assert two_equal_contigs != unequal_contigs


def test_prokanota_id_preserves_sequence_order():
    """
    Verify reordering otherwise identical contigs changes the ID.
    Prokanota assigns feature numbers in input order, so sequence order is
    part of the identifier and must be represented by the hash.
    """
    original_order = feature_utils.generate_prokanota_id(
        sample_id="test_sample",
        mode="genome",
        sequences=["ACGT", "TGCA"],
    )
    reversed_order = feature_utils.generate_prokanota_id(
        sample_id="test_sample",
        mode="genome",
        sequences=["TGCA", "ACGT"],
    )

    assert original_order != reversed_order


def test_prokanota_id_separates_genome_and_protein_modes():
    """
    Verify identical character sequences differ between input modes.
    A sequence such as ACGT is valid text for both DNA and protein input, so
    the mode must prevent cross-mode identifier reuse.
    """
    genome_id = feature_utils.generate_prokanota_id(
        sample_id="test_sample",
        mode="genome",
        sequences=["ACGT"],
    )
    protein_id = feature_utils.generate_prokanota_id(
        sample_id="test_sample",
        mode="protein",
        sequences=["ACGT"],
    )

    assert genome_id != protein_id


def test_prokanota_id_normalizes_sequence_case():
    """
    Verify upper- and lowercase FASTA sequences produce the same ID.
    Sequence case is formatting rather than a biological difference.
    """
    uppercase_id = feature_utils.generate_prokanota_id(
        sample_id="test_sample",
        mode="genome",
        sequences=["ACGT"],
    )
    lowercase_id = feature_utils.generate_prokanota_id(
        sample_id="test_sample",
        mode="genome",
        sequences=["acgt"],
    )

    assert uppercase_id == lowercase_id


def test_prokanota_id_ignores_terminal_protein_stop_symbol():
    """
    Verify an optional terminal protein stop symbol does not change the ID.
    The terminal star is translation notation rather than an amino acid and
    is removed from Prokanota's imported protein output.
    """
    without_stop = feature_utils.generate_prokanota_id(
        sample_id="test_sample",
        mode="protein",
        sequences=["MPEPTIDE"],
    )
    with_stop = feature_utils.generate_prokanota_id(
        sample_id="test_sample",
        mode="protein",
        sequences=["MPEPTIDE*"],
    )

    assert without_stop == with_stop


def test_prokanota_id_rejects_unknown_mode():
    """
    Verify misspelled or unsupported input modes are rejected.
    Accepting an unknown mode could create IDs outside the documented genome
    and protein namespaces.
    """
    with pytest.raises(ValueError, match="genome.*protein"):
        feature_utils.generate_prokanota_id(
            sample_id="test_sample",
            mode="unknown",
            sequences=["ACGT"],
        )


# ---
# 7. Sample ID Validation (validate_metadata_csv)
# ---


def test_duplicate_sample_ids_are_rejected(tmp_path):
    """Duplicate sample IDs in metadata should raise an error."""
    pytest.importorskip("snaketool_utils")
    import click

    from prokanota.__main__ import validate_metadata_csv

    metadata_path = tmp_path / "metadata.csv"
    metadata_path.write_text(
        "sample_id,path\nsample_1,first.fasta\nsample_1,second.fasta\n",
        encoding="utf-8",
    )

    with pytest.raises(click.ClickException, match="duplicate sample_id"):
        validate_metadata_csv(metadata_path)
