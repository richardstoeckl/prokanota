"""Biologically grounded tests for published truths."""

from __future__ import annotations

from prokanota.workflow.scripts import feature_utils

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
    UniProt lists MW as 90390 Da, ExPASy ProtParam calculates 90390.15 Da.
    """
    pf_ago = "MKAKVVINLVKINKKIIPDKIYVYRLFNDPEEELQKEGYSIYRLAYENVGIVIDPENLIIATTKELEYEGEFIPEGEISFSELRNDYQSKLVLRLLKENGIGEYELSKLLRKFRKPKTFGDYKVIPSVEMSVIKHDEDFYLVIHIIHQIQSMKTLWELVNKDPKELEEFLMTHKENLMLKDIASPLKTVYKPCFEEYTKKPKLDHNQEIVKYWYNYHIERYWNTPEAKLEFYRKFGQVDLKQPAILAKFASKIKKNKNYKIYLLPQLVVPTYNAEQLESDVAKEILEYTKLMPEERKELLENILAEVDSDIIDKSLSEIEVEKIAQELENKIRVRDDKGNSVPISQLNVQKSQLLLWTNYSRKYPVILPYEVPEKFRKIREIPMFIILDSGLLADIQNFATNEFRELVKSMYYSLAKKYNSLAKKARSTNEIGLPFLDFRGKEKVITEDLNSDKGIIEVVEQVSSFMKGKELGLAFIAARNKLSSEKFEEIKRRLFNLNVISQVVNEDTLKNKRDKYDRNRLDLFVRHNLLFQVLSKLGVKYYVLDYRFNYDYIIGIDVAPMKRSEGYIGGSAVMFDSQGYIRKIVPIKIGEQRGESVDMNEFFKEMVDKFKEFNIKLDNKKILLLRDGRITNNEEEGLKYISEMFDIEVVTMDVIKNHPVRAFANMKMYFNLGGAIYLIPHKLKQAKGTPIPIKLAKKRIIKNGKVEKQSITRQDVLDIFILTRLNYGSISADMRLPAPVHYAHKFANAIRNEWKIKEEFLAEGFLYFV"
    mw = feature_utils.protein_molecular_weight(pf_ago)
    assert abs(mw - 90390) < 1.0  # allow for small rounding differences


# ---
# 2. DNA Reverse Complement (feature_utils.py)
# ---

# Reverse complementation is fundamental for reporting features predicted on
# the minus strand. Canonical bases follow Watson-Crick pairing (A<->T, C<->G),
# while ambiguous bases follow the IUPAC nucleotide-code definitions:
# https://www.bioinformatics.org/sms/iupac.html


def test_reverse_complement_basic():
    """
    Verify basic reverse complement transformation.
    Input:  5'-ATGC-3'
    Output: 5'-GCAT-3'
    """
    assert feature_utils.reverse_complement("ATGC") == "GCAT"


def test_reverse_complement_palindrome():
    """
    Test a palindromic sequence using the EcoRI recognition site.
    GAATTC equals its own reverse complement.
    """
    ecori_site = "GAATTC"
    assert feature_utils.reverse_complement(ecori_site) == ecori_site


def test_reverse_complement_iupac_bases():
    """
    Verify complementation of IUPAC ambiguity codes.

    R represents A or G and therefore complements to Y (C or T). Because the
    complete sequence is also reversed, 5'-ATRC-3' becomes 5'-GYAT-3'.
    """
    assert feature_utils.reverse_complement("ATRC") == "GYAT"


def test_single_nucleotide_reverse_complement():
    """Verify the canonical Watson-Crick complement of each DNA nucleotide."""
    assert feature_utils.reverse_complement("A") == "T"
    assert feature_utils.reverse_complement("T") == "A"
    assert feature_utils.reverse_complement("G") == "C"
    assert feature_utils.reverse_complement("C") == "G"
