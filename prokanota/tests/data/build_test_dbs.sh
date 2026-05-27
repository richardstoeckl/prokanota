#!/usr/bin/env bash
# =============================================================================
# build_test_dbs.sh
# Build HMM, BLAST, DIAMOND, MMseqs2, and RPS-BLAST test databases from
# arCOG alignments, and produce annotation mapping files.
#
# Dependencies: wget, tar, hmmbuild, hmmpress, esl-reformat (HMMER3),
#               seqkit, makeblastdb, psiblast, makeprofiledb, diamond, mmseqs
#
# =============================================================================
# This script is explicitly NOT intended to be used by prokanota itself; it is 
# a utility for building the test databases used in prokanota's own test suite.
# The resulting databases are stored in the repository under tests/data/mockDB/
# and are not designed to be used outside of the test suite. They are built
# from an arbitrary small subset of arCOG profiles and proteins. The script is
# included in the repository for transparency and reproducibility only.
# =============================================================================
#
# Output layout:
#   logs/               — all log files
#   interim/            — downloaded archive, .sr files, intermediate files
#   mockDB/
#     hmm/              — hmm_mockDB.hmm + .h3f/.h3i/.h3m/.h3p
#     blastp/           — blastp_mockDB.p*
#     diamond/          — diamond_mockDB.dmnd
#     mmseqs2/          — mmseqs2_mockDB (MMseqs2 index files)
#     rpsblast/         — rpsblast_mockDB.rps/.loo/.aux/...
#     mockDB_profile_mapping.tsv — per-profile annotation (arCOG ID → name/desc/cat)
#     mockDB_seq_mapping.tsv     — per-sequence annotation (seq ID → name/desc/cat)
#
# All sequence databases (blastp, diamond, mmseqs2, rpsblast) are built from
# the exact same source: interim/all_seqs.faa (degapped, cleaned protein FASTA).
#
# Usage: bash build_test_dbs.sh
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
ARCOG_IDS=(
    arCOG03890 arCOG04089 arCOG01452 arCOG03422 arCOG01829 arCOG01920
)

ARCOG_URL="https://ftp.ncbi.nlm.nih.gov/pub/wolf/COGs/arCOG/zip.aliar14.tgz"
ARCOGDEF_URL="https://ftp.ncbi.nlm.nih.gov/pub/wolf/COGs/arCOG/ar14.arCOGdef.tab"

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------
LOGS_DIR="logs"
INTERIM_DIR="interim"
MOCKDB_DIR="mockDB"

SR_DIR="${INTERIM_DIR}/sr"           # extracted .sr alignment files
AFA_DIR="${INTERIM_DIR}/afa"         # cleaned aligned FASTA files
SMP_DIR="${INTERIM_DIR}/smp"         # PSSM scoremat files (for RPS step)
HMM_PARTS_DIR="${INTERIM_DIR}/hmm"   # per-profile .hmm files before concatenation

HMM_OUT="${MOCKDB_DIR}/hmm"
BLASTP_OUT="${MOCKDB_DIR}/blastp"
DIAMOND_OUT="${MOCKDB_DIR}/diamond"
MMSEQS2_OUT="${MOCKDB_DIR}/mmseqs2"
RPS_OUT="${MOCKDB_DIR}/rpsblast"

HMM_DB="${HMM_OUT}/hmm_mockDB"
BLASTP_DB="${BLASTP_OUT}/blastp_mockDB"
DIAMOND_DB="${DIAMOND_OUT}/diamond_mockDB"
MMSEQS2_DB="${MMSEQS2_OUT}/mmseqs2_mockDB"
RPS_DB="${RPS_OUT}/rpsblast_mockDB"

ARCHIVE="${INTERIM_DIR}/zip.aliar14.tgz"
ARCOGDEF_RAW="${INTERIM_DIR}/arCOGdef.tab"
SEQS_FASTA="${INTERIM_DIR}/all_seqs.faa"   # single source of truth for all
                                            # sequence databases — degapped,
                                            # cleaned protein FASTA
SMP_LIST="${INTERIM_DIR}/smp_list.pn"
PROFILE_MAPPING="${MOCKDB_DIR}/mockDB_profile_mapping.tsv"
SEQ_MAPPING="${MOCKDB_DIR}/mockDB_seq_mapping.tsv"

# ---------------------------------------------------------------------------
# HELPER: check required tools
# ---------------------------------------------------------------------------
check_deps() {
    local missing=()
    for cmd in wget tar hmmbuild hmmpress esl-reformat seqkit \
               makeblastdb psiblast makeprofiledb diamond mmseqs; do
        command -v "$cmd" &>/dev/null || missing+=("$cmd")
    done
    if [[ ${#missing[@]} -gt 0 ]]; then
        echo "ERROR: missing dependencies: ${missing[*]}" >&2
        exit 1
    fi
}

check_deps
mkdir -p "${LOGS_DIR}" "${SR_DIR}" "${AFA_DIR}" "${SMP_DIR}" \
         "${HMM_PARTS_DIR}" "${HMM_OUT}" "${BLASTP_OUT}" \
         "${DIAMOND_OUT}" "${MMSEQS2_OUT}" "${RPS_OUT}"

# ===========================================================================
# STEP 1 — Download and extract selected .sr alignment files
#
# Files are stored as plain (uncompressed) .sr files under "ar14.ali/" inside
# the archive. All paths are passed to a single tar call so the archive is
# only scanned once. --strip-components=1 drops the "ar14.ali/" prefix.
# ===========================================================================
if [[ ! -f "${ARCHIVE}" ]]; then
    echo "[1/6] Downloading arCOG alignment archive..."
    wget -q --show-progress -O "${ARCHIVE}" "${ARCOG_URL}"
else
    echo "[1/6] Archive already present, skipping download."
fi

echo "[1/6] Extracting ${#ARCOG_IDS[@]} selected .sr files..."

TAR_PATHS=()
for id in "${ARCOG_IDS[@]}"; do
    TAR_PATHS+=( "ar14.ali/${id}.sr" )
done

tar -xzf "${ARCHIVE}" \
    -C "${SR_DIR}" \
    --strip-components=1 \
    "${TAR_PATHS[@]}"

SR_FILES=( "${SR_DIR}/"*.sr )
echo "    Extracted ${#SR_FILES[@]} .sr file(s)."

# ===========================================================================
# STEP 2 — Build a pressed HMM profile database (for pyhmmer)
#
# hmmbuild reads SELEX (.sr) format natively and builds one profile per file.
# All profiles are concatenated, then hmmpress creates the binary index files
# (.h3f/.h3i/.h3m/.h3p) that pyhmmer requires for optimised search.
# ===========================================================================
echo "[2/6] Building HMM database..."

HMM_CAT="${HMM_DB}.hmm"
> "${HMM_CAT}"

for sr in "${SR_FILES[@]}"; do
    id=$(basename "${sr}" .sr)
    hmmbuild --amino -n "${id}" \
        "${HMM_PARTS_DIR}/${id}.hmm" "${sr}" \
        > "${LOGS_DIR}/${id}.hmmbuild.log" 2>&1
    cat "${HMM_PARTS_DIR}/${id}.hmm" >> "${HMM_CAT}"
done

hmmpress "${HMM_CAT}" > "${LOGS_DIR}/hmmpress.log" 2>&1
echo "    HMM database: ${HMM_DB}.hmm (+ .h3f/.h3i/.h3m/.h3p)"

# ===========================================================================
# STEP 3 — Convert alignments to clean aligned FASTA and build the shared
#           protein FASTA (interim/all_seqs.faa)
#
# This file is the single source of truth for all sequence-level databases
# (blastp, diamond, mmseqs2, rpsblast). Building them all from the same FASTA
# guarantees identical sequence content across tools.
#
# psiblast -in_msa's CAlnReader is strict about its aligned FASTA input.
# Three classes of input defect all produce the same misleading error:
#   "CAlnReader::GetSeqEntry(): Seq_entry is not available until after Read()"
#
# Defect 1 — Non-alphanumeric characters in sequence NAMES (e.g. '|', ';'):
#   Fix: use '_' as the separator when prepending the arCOG ID.
#
# Defect 2 — Non-standard residue characters in SEQUENCES ('*', '#', '+', '.'):
#   arCOG alignments use '*' for stop codons and '.' for insert-column gaps.
#   Fix: map all four to '-' with tr. Every sequence is preserved intact;
#   only the representation of those positions changes to a standard gap.
#
# Defect 3 — Duplicate sequence NAMES within a single MSA:
#   arCOG alignments can contain paralogous sequences with identical IDs.
#   Fix: seqkit rename (NOT seqkit rmdup). seqkit rename appends _N to make
#   each name unique while keeping every sequence in the file. seqkit rmdup
#   would silently discard sequences and must never be used here.
# ===========================================================================
echo "[3/6] Cleaning alignments and building shared protein FASTA..."

> "${SEQS_FASTA}"

for sr in "${SR_FILES[@]}"; do
    id=$(basename "${sr}" .sr)

    # Convert SELEX → aligned FASTA, then:
    #   seqkit replace  : prepend "<arCOG_id>_" (underscore avoids CAlnReader
    #                     rejection of special chars in names)
    #   seqkit rename   : make duplicate names unique by appending _1/_2/...
    #                     (keeps every sequence — never drops anything)
    #   seqkit seq -u   : uppercase residues
    #   tr '*#+.' '----': map stop-codon / insert-gap symbols to standard gaps
    esl-reformat afa "${sr}" \
        | seqkit replace -p "^" -r "${id}_" \
        | seqkit rename \
        | seqkit seq -u \
        | tr '*#+.' '----' \
        > "${AFA_DIR}/${id}.afa"

    # Degap → plain protein sequences, appended to the shared FASTA
    seqkit seq -g "${AFA_DIR}/${id}.afa" >> "${SEQS_FASTA}"
done

NSEQS=$(grep -c '^>' "${SEQS_FASTA}")
echo "    Shared protein FASTA: ${SEQS_FASTA} (${NSEQS} sequences)"

# ===========================================================================
# STEP 4 — Build sequence-level databases from the shared FASTA
#
# All three databases are built from the exact same ${SEQS_FASTA}, so any
# search result differences between tools are due to the tools themselves,
# not to differences in the underlying sequences.
#
# 4a. BLAST protein database
#     -parse_seqids is required so psiblast (step 5) can look up sequences
#     by name when computing PSSMs.
#
# 4b. DIAMOND database
#     diamond makedb --in <fasta> --db <output>
#     Produces a single binary .dmnd file.
#
# 4c. MMseqs2 database
#     mmseqs createdb <fasta> <output_prefix>
#     Produces a set of flat binary index files (no conversion from BLAST
#     format is possible — MMseqs2 always requires its own format).
# ===========================================================================
echo "[4/6] Building BLAST, DIAMOND, and MMseqs2 databases..."

# 4a — BLAST
makeblastdb \
    -in "${SEQS_FASTA}" \
    -dbtype prot \
    -out "${BLASTP_DB}" \
    -parse_seqids \
    -title "blastp_mockDB" \
    > "${LOGS_DIR}/makeblastdb.log" 2>&1
echo "    BLAST database:   ${BLASTP_DB}.p*"

# 4b — DIAMOND
diamond makedb \
    --in "${SEQS_FASTA}" \
    --db "${DIAMOND_DB}" \
    > "${LOGS_DIR}/diamond_makedb.log" 2>&1
echo "    DIAMOND database: ${DIAMOND_DB}.dmnd"

# 4c — MMseqs2
mmseqs createdb \
    "${SEQS_FASTA}" \
    "${MMSEQS2_DB}" \
    > "${LOGS_DIR}/mmseqs2_createdb.log" 2>&1
echo "    MMseqs2 database: ${MMSEQS2_DB}"

# ===========================================================================
# STEP 5 — Build an RPS-BLAST database
#
# Uses the cleaned aligned-FASTA files from step 3 (already free of the
# characters that crash CAlnReader, and with arCOG-prefixed sequence names
# matching the entries in the BLAST db built in step 4a).
#
# psiblast -in_msa reads each .afa and computes a PSSM scoremat (.smp) from
# the column frequencies of the MSA. The -db argument is required by psiblast
# for sequence lookup; we reuse the blastp database from step 4a.
# The search output is discarded (-out /dev/null); only -out_pssm matters.
#
# makeprofiledb assembles all .smp files into the final RPS-BLAST database.
# ===========================================================================
echo "[5/6] Building RPS-BLAST database..."

> "${SMP_LIST}"

for sr in "${SR_FILES[@]}"; do
    id=$(basename "${sr}" .sr)
    smp="${SMP_DIR}/${id}.smp"

    psiblast \
        -in_msa "${AFA_DIR}/${id}.afa" \
        -db "${BLASTP_DB}" \
        -num_iterations 1 \
        -out_pssm "${smp}" \
        -out /dev/null \
        > "${LOGS_DIR}/${id}.psiblast.log" 2>&1

    echo "${smp}" >> "${SMP_LIST}"
done

makeprofiledb \
    -in "${SMP_LIST}" \
    -out "${RPS_DB}" \
    -title "rpsblast_mockDB" \
    -dbtype rps \
    > "${LOGS_DIR}/makeprofiledb.log" 2>&1

echo "    RPS-BLAST database: ${RPS_DB}.rps/.loo/.aux/..."

# ===========================================================================
# STEP 6 — Download arCOG definition file and produce mapping tables
#
# Source: ar14.arCOGdef.tab
# Column layout (tab-separated, 1-indexed):
#   $1  arCOG ID
#   $2  functional category letter(s)
#   $3  arCOG name / short description  (may contain '-' padding artefacts)
#   $4  longer functional description
#
# Two output files are produced:
#
# 6a. mockDB_profile_mapping.tsv  — one row per arCOG profile
#       columns: arCOG_ID | cleaned_name | description | category
#       Used directly for HMM search result annotation (query name = arCOG ID).
#
# 6b. mockDB_seq_mapping.tsv  — one row per individual protein sequence
#       columns: seq_ID | cleaned_name | description | category
#       seq_ID has the form "<arCOG_ID>_<accession>" (e.g. arCOG01829_80000072),
#       exactly matching the sequence names in the BLAST/DIAMOND/MMseqs2 databases.
#       Built by extracting all headers from all_seqs.faa, then joining on the
#       arCOG prefix (everything up to and including the last '_<accession>'
#       suffix, i.e. the part before the first '_8...' numeric accession).
#
#       Join logic (pure awk, no temp files, no extra tools):
#         Pass 1 (-v ids=...): read arCOGdef.tab, build profile_map[arCOG_ID]
#                              = "name\tdesc\tcat" for the selected IDs only.
#         Pass 2:              read FASTA headers (lines starting with '>'),
#                              strip '>', extract the arCOG prefix by removing
#                              the trailing '_<numeric_accession>' suffix using
#                              a greedy match on the last underscore-delimited
#                              numeric field, then look up in profile_map.
#
#       The arCOG prefix extraction uses the pattern:
#           sub(/_[0-9]+(_[0-9]+)*$/, "", prefix)
#       which strips the trailing numeric accession added by seqkit rename
#       (e.g. _80000072, or _80000072_1 for renamed duplicates), leaving the
#       bare arCOG ID (e.g. arCOG01829).
# ===========================================================================
echo "[6/6] Downloading arCOG definition file and building mapping tables..."

if [[ ! -f "${ARCOGDEF_RAW}" ]]; then
    wget -q --show-progress -O "${ARCOGDEF_RAW}" "${ARCOGDEF_URL}"
else
    echo "    Definition file already present, skipping download."
fi

# Pipe-separated ID string for awk filtering (avoids shell quoting issues)
IDS_PIPE=$(IFS='|'; echo "${ARCOG_IDS[*]}")

# ---------------------------------------------------------------------------
# 6a — Profile mapping (one row per arCOG, identical to previous version)
# ---------------------------------------------------------------------------
tail -n +2 "${ARCOGDEF_RAW}" \
    | LC_CTYPE=C sed 's/\r//g' \
    | awk -F'\t' -v ids="${IDS_PIPE}" \
        'BEGIN {
            OFS = "\t"
            n = split(ids, a, "|")
            for (i = 1; i <= n; i++) keep[a[i]] = 1
        }
        {
            if (!($1 in keep)) next
            gsub(/-/, "", $3)
            print $1, $3, $4, $2
        }' \
    > "${PROFILE_MAPPING}"

NPROFILES=$(wc -l < "${PROFILE_MAPPING}")
echo "    Profile mapping: ${PROFILE_MAPPING} (${NPROFILES} entries)"

# ---------------------------------------------------------------------------
# 6b — Sequence mapping (one row per protein sequence in the databases)
#
# Two-pass awk over two input files:
#   File 1: arCOGdef.tab  → build profile_map[arCOG_ID] = "name\tdesc\tcat"
#   File 2: all_seqs.faa  → for each '>' header, derive arCOG prefix and look up
# ---------------------------------------------------------------------------
awk -F'\t' -v ids="${IDS_PIPE}" \
    'BEGIN {
        OFS = "\t"
        # Build lookup set of wanted IDs
        n = split(ids, a, "|")
        for (i = 1; i <= n; i++) keep[a[i]] = 1
    }

    # --- Pass 1: arCOGdef.tab ---
    # NR == FNR is true only while reading the first file
    NR == FNR {
        if (FNR == 1) next          # skip header
        id = $1
        if (!(id in keep)) next
        name = $3; gsub(/-/, "", name)
        profile_map[id] = name "\t" $4 "\t" $2
        next
    }

    # --- Pass 2: all_seqs.faa ---
    # Only process FASTA header lines
    /^>/ {
        seq_id = substr($0, 2)      # strip leading ">"
        # Derive the arCOG prefix by stripping the trailing numeric accession(s).
        # Sequence names have the form:  arCOG01829_80000072
        # or after seqkit rename:        arCOG01829_80000072_1
        # Strategy: remove the last underscore-separated token if it is purely
        # numeric, repeating until the remaining string is a known arCOG ID.
        # A single sub() on "_[0-9]+$" handles both cases because seqkit rename
        # appends _1/_2/... as a purely numeric suffix after the accession.
        prefix = seq_id
        sub(/_[0-9]+$/, "", prefix)     # strip rename suffix (_1, _2, ...) if present
        sub(/_[0-9]+$/, "", prefix)     # strip numeric accession (e.g. _80000072)
        if (prefix in profile_map) {
            print seq_id, profile_map[prefix]
        }
    }' \
    "${ARCOGDEF_RAW}" "${SEQS_FASTA}" \
    > "${SEQ_MAPPING}"

NSEQ_ROWS=$(wc -l < "${SEQ_MAPPING}")
echo "    Sequence mapping: ${SEQ_MAPPING} (${NSEQ_ROWS} entries)"

# ===========================================================================
# DONE
# ===========================================================================
echo ""
echo "All databases and mapping tables built successfully:"
echo ""
echo "  mockDB/"
echo "    hmm/                       ${HMM_DB}.hmm  (+ .h3f .h3i .h3m .h3p)"
echo "    blastp/                    ${BLASTP_DB}.p*"
echo "    diamond/                   ${DIAMOND_DB}.dmnd"
echo "    mmseqs2/                   ${MMSEQS2_DB}"
echo "    rpsblast/                  ${RPS_DB}.rps .loo .aux ..."
echo "    mockDB_profile_mapping.tsv ${NPROFILES} arCOG profiles"
echo "    mockDB_seq_mapping.tsv     ${NSEQ_ROWS} protein sequences"
echo ""
echo "  Shared sequence source: ${SEQS_FASTA} (${NSEQS} sequences)"
echo "  logs/    $(ls "${LOGS_DIR}"/ | wc -l) log files"
echo "  interim/ source .sr, .afa, .smp files, individual .hmm files"
