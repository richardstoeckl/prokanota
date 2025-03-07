"""

Copyright Richard Stöckl 2025.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE or copy at 
https://www.boost.org/LICENSE_1_0.txt)

"""

import hashlib
import os
import pyrodigal
import argparse
import logging
import gb_io
from gb_io import Record, Feature, Qualifier

# Configure logging at the INFO level
logging.basicConfig(level=logging.INFO)
log = logging.getLogger('CDS')

class GenePrediction:
    """
    A container class for storing information about a gene prediction.

    Attributes:
        contig_id (str): Identifier for the contig the gene is on.
        gene_id (str): Unique ID of the gene (includes a zero-padded index).
        start (int): Start coordinate of the gene on the contig.
        end (int): End coordinate of the gene on the contig.
        strand (str): '+' or '-' denoting the gene's strand.
        protein_seq (str): The amino acid sequence of the predicted gene.
        dna_seq (str): The nucleotide sequence of the predicted gene.
    """
    def __init__(self, contig_id, gene_id, start, end, strand, protein_seq, dna_seq):
        self.contig_id = contig_id
        self.gene_id = gene_id
        self.start = start
        self.end = end
        self.strand = strand
        self.protein_seq = protein_seq
        self.dna_seq = dna_seq

def parse_fasta_and_predict(sample_id, filename, meta_mode, closed, translation_table):
    """
    Reads the FASTA file, computes a genome ID using MD5 hashing and a custom mapping,
    then uses Pyrodigal to predict genes.
    
    Additionally, creates a dictionary mapping contig_tag to the original FASTA header id.
    
    Returns:
        tuple: (genome_id, tagged_contigs, gene_records, contig_mapping)
            - genome_id (str): The generated 8-character genome ID.
            - tagged_contigs (dict): A dictionary of contig_tag -> contig_sequence.
            - gene_records (list): A list of GenePrediction objects.
            - contig_mapping (dict): Mapping of contig_tag -> original FASTA header id.
    """
    # Initialize MD5 hash object with the sample_id
    hash_obj = hashlib.md5(sample_id.encode('utf-8'))

    # Define the mapping from hex digits (0-9,a-f) to letters (A-J,K-P)
    mapping = {str(i): chr(65 + i) for i in range(10)}
    mapping.update({chr(97 + i): chr(75 + i) for i in range(6)})

    # Read the FASTA file to build a dict {original_header: sequence}, updating the hash
    contigs = {}
    current_header, current_sequence = None, []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_header is not None:
                    contigs[current_header] = ''.join(current_sequence)
                current_header = line.strip()[1:]
                current_sequence = []
            else:
                seq_line = line.strip().upper()
                hash_obj.update(seq_line.encode('utf-8'))
                current_sequence.append(seq_line)
        if current_header is not None:
            contigs[current_header] = ''.join(current_sequence)

    # Compute the genome ID from the first 8 hex digits
    hexdigest = hash_obj.hexdigest()
    genome_id = ''.join(mapping[c] for c in hexdigest[:8])

    # Create new dictionaries: one mapping contig_tag -> sequence and one for contig_tag -> original header id
    tagged_contigs = {}
    contig_mapping = {}
    for i, (header, seq) in enumerate(contigs.items()):
        contig_tag = f"{genome_id}_{i+1}"
        # Extrahiere nur die erste Spalte des Headers (ohne Beschreibung)
        header_id = header.split()[0]
        tagged_contigs[contig_tag] = seq
        contig_mapping[contig_tag] = header_id

    # Create a Pyrodigal GeneFinder. Train if meta_mode is off
    gene_finder = pyrodigal.GeneFinder(meta=meta_mode, closed=closed)
    if not meta_mode:
        gene_finder.train(*tagged_contigs.values(), translation_table=translation_table)

    gene_records = []
    # Iterate through each tagged contig and predict genes
    for contig_tag, dna_sequence in tagged_contigs.items():
        sequence = pyrodigal.Sequence(dna_sequence)
        genes = gene_finder.find_genes(sequence)
        for j, gene in enumerate(genes):
            protein_seq = gene.translate()
            gene_id = f"{contig_tag}_{j+1:05d}"
            strand_symbol = '+' if gene.strand == 1 else '-'
            gene_dna_seq = dna_sequence[gene.begin-1:gene.end]
            gene_records.append(
                GenePrediction(
                    contig_id=contig_tag,
                    gene_id=gene_id,
                    start=gene.begin,
                    end=gene.end,
                    strand=strand_symbol,
                    protein_seq=protein_seq,
                    dna_seq=gene_dna_seq
                )
            )

    return genome_id, tagged_contigs, gene_records, contig_mapping

def write_faa(genome_id, gene_records, faa_path):
    """
    Writes all predicted proteins to a FASTA (.faa) file.
    Each record is written as a header containing the gene_id,
    followed by the amino acid sequence.

    Args:
        genome_id (str): The unique genome ID, used as the file name prefix.
        gene_records (list): A list of GenePrediction objects.
        faa_path (str): Path to the output .faa file.
    """
    with open(faa_path, "w") as faa_file:
        for record in gene_records:
            faa_file.write(f">{record.gene_id}\n{record.protein_seq}\n")

def write_gff(genome_id, gene_records, gff_path):
    """
    Writes all predicted gene locations to a GFF file.
    Each record includes the contig_id, start, end, and gene ID as an attribute.

    Args:
        genome_id (str): The unique genome ID, used as the file name prefix.
        gene_records (list): A list of GenePrediction objects.
        gff_path (str): Path to the output .gff file.
    """
    with open(gff_path, "w") as gff_file:
        for record in gene_records:
            gff_file.write(
                f"{record.contig_id}\tPyrodigal\tgene\t"
                f"{record.start}\t{record.end}\t.\t"
                f"{record.strand}\t.\tID={record.gene_id}\n"
            )

def write_gbk(genome_id, contigs, gene_records, contig_mapping, gbk_path):
    """
    Writes the predicted CDS to a GenBank file.

    Args:
        genome_id (str): The generated 8-character genome ID.
        contigs (dict): A dictionary of contig_header -> contig_sequence.
        gene_records (list): A list of GenePrediction objects.
        contig_mapping (dict): Mapping of contig_tag -> original FASTA header id.
        gbk_path (str): The output path for the GenBank file.
    """
    records = []
    for contig_tag, header_id in contig_mapping.items():
        # Hole die Sequenz über contig_tag
        dna_sequence = contigs[contig_tag]
        # Encode to bytes if required
        if isinstance(dna_sequence, str):
            dna_sequence = dna_sequence.encode('utf-8')

        features = []
        for gene in gene_records:
            if gene.contig_id == contig_tag:
                qualifiers = [
                    Qualifier(key="locus_tag", value=gene.gene_id),
                    Qualifier(key="translation", value=gene.protein_seq)
                ]
                # Adjust the start for plus strand features only.
                if gene.strand == '+':
                    adjusted_start = gene.start - 1
                else:
                    adjusted_start = gene.start

                base_location = gb_io.Range(start=adjusted_start, end=gene.end)
                # For negative strand, wrap in Complement.
                if gene.strand == '-':
                    location = gb_io.Complement(base_location)
                else:
                    location = base_location

                feature = Feature(
                    kind="CDS",
                    location=location,
                    qualifiers=qualifiers
                )
                features.append(feature)


        record = Record(
            name=header_id,
            sequence=dna_sequence,
            features=features,
            accession=genome_id
        )
        records.append(record)

    with open(gbk_path, "wb") as file:
        gb_io.dump(records, file)


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]

def write_fna(genome_id, gene_records, fna_path):
    """
    Writes gene sequences to a FASTA (.fna) file.
    For genes on the plus strand, uses the sequence as predicted.
    For genes on the minus strand, reverse-complements the predicted sequence,
    then trims the last nucleotide to correct for coordinate discrepancies.
    """
    with open(fna_path, "w") as fna_file:
        for record in gene_records:
            if record.strand == '-':
                # Reverse complement the entire gene sequence,
                # then trim the last nucleotide.
                corrected_seq = reverse_complement(record.dna_seq)[:-1]
            else:
                corrected_seq = record.dna_seq
            fna_file.write(f">{record.gene_id}\n{corrected_seq}\n")


def protein_molecular_weight(seq):
    """
    Estimate the molecular weight of a protein sequence in Daltons.
    Uses average residue masses (in Da) and adds the mass of one water molecule.
    This method is described in GASTEIGER, Elisabeth, et al. The proteomics protocols handbook, 2005, S. 571-607.
    The values can be found at https://web.expasy.org/findmod/findmod_masses.html#AA
    """
    weights = {
      'A': 71.0788, 'R': 156.1875, 'N': 114.1038, 'D': 115.0886, 'C': 103.1388,
      'E': 129.1155, 'Q': 128.1307, 'G': 57.0513, 'H': 137.1411, 'I': 113.1594,
      'L': 113.1594, 'K': 128.1741, 'M': 131.1926, 'F': 147.1766, 'P': 97.1167,
      'S': 87.0782, 'T': 101.1051, 'W': 186.2132, 'Y': 163.1733, 'V': 99.1311
    }
    mw = 0.0
    for aa in seq:
        mw += weights.get(aa, 0.0)  # if an unknown amino acid appears, weight is 0 to alert the user
    # Add the weight of one water molecule (approximately 18.015 Da) for the free termini
    mw += 18.01528
    return mw



def write_tsv(sample_id, contigs, gene_records, contig_mapping, tsv_path):
    """
    Writes the gene predictions to a TSV file.

    Args:
        sample_id (str): The sample ID.
        contigs (dict): A dictionary of contig_header -> contig_sequence.
        gene_records (list): A list of GenePrediction objects.
        contig_mapping (dict): Mapping of contig_tag -> original FASTA header id.
        tsv_path (str): The output path for the TSV file.
    """
    with open(tsv_path, "w") as tsv_file:
        # Write header with the new columns
        tsv_file.write("sample_id\torig_cont_header\tcontig_id\tcontig_length\t"
                       "gene_id\tstart\tend\tstrand\tgene_length\tprotein_length\tprotein_mw_kDa\n")
        for contig_tag, header_id in contig_mapping.items():
            contig_length = len(contigs[contig_tag])
            for gene in gene_records:
                if gene.contig_id == contig_tag:
                    # Calculate gene length. For genes on '-' strand, subtract one nucleotide (see write_gbk()).
                    gene_length = len(gene.dna_seq)
                    if gene.strand == '-':
                        gene_length -= 1

                    # Remove stop codon from the protein sequence if present for length calculation
                    if gene.protein_seq.endswith("*"):
                        protein_seq = gene.protein_seq[:-1]
                    else:
                        protein_seq = gene.protein_seq
                    protein_length = len(protein_seq)

                    # Calculate molecular weight (in Daltons),
                    # then convert to kDa and format to one decimal place for output
                    mw_da = protein_molecular_weight(protein_seq)
                    mw_kda = mw_da / 1000.0

                    tsv_file.write(
                        f"{sample_id}\t{header_id}\t{contig_tag}\t{contig_length}\t"
                        f"{gene.gene_id}\t{gene.start}\t{gene.end}\t{gene.strand}\t"
                        f"{gene_length}\t{protein_length}\t{mw_kda:.1f}\n"
                    )
    log.info(f"TSV summary written to {tsv_path}")

def process_genome(sample_id, filename, faa_path, gff_path, gbk_path, fna_path, tsv_path, 
                   meta_mode, closed, translation_table, 
                   write_faa_flag, write_gff_flag, write_gbk_flag, write_fna_flag, write_tsv_flag):
    """
    Master function to process a single genome. Parses the FASTA, runs predictions,
    and writes .faa, .gff, .gbk, .fna files, and .tsv files based on flags.

    Args:
        sample_id (str): The user-provided sample ID.
        filename (str): Path to the genome FASTA file.
        faa_path (str): Output path for the predicted .faa file.
        gff_path (str): Output path for the .gff file.
        gbk_path (str): Output path for the .gbk file.
        fna_path (str): Output path for the .fna file.
        tsv_path (str): Output path for the .tsv file.
        meta_mode (bool): True if running Pyrodigal in metagenomic mode.
        closed (bool): True if the genome is closed-ended.
        translation_table (int): Translation table used to decode the DNA into proteins.
        write_faa_flag (bool): Flag to write .faa file.
        write_gff_flag (bool): Flag to write .gff file.
        write_gbk_flag (bool): Flag to write .gbk file.
        write_fna_flag (bool): Flag to write .fna file.
        write_tsv_flag (bool): Flag to write .tsv file.
    """
    # Parse the FASTA and get predictions
    genome_id, contigs, gene_records, contig_mapping = parse_fasta_and_predict(
        sample_id, filename, meta_mode, closed, translation_table
    )

    # Write outputs based on flags
    if write_faa_flag:
        write_faa(genome_id, gene_records, faa_path)
    if write_gff_flag:
        write_gff(genome_id, gene_records, gff_path)
    if write_gbk_flag:
        write_gbk(genome_id, contigs, gene_records, contig_mapping, gbk_path)
    if write_fna_flag:
        write_fna(genome_id, gene_records, fna_path)
    if write_tsv_flag:
        write_tsv(sample_id, contigs, gene_records, contig_mapping, tsv_path)

    # Log processed file info
    log.info(f"Processed {filename} => {genome_id}")

if __name__ == "__main__":
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(
        description="Predict genes and write results to .faa, .gff, .gbk, .fna, and .tsv files."
    )
    parser.add_argument("sample_id", help="Sample ID")
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument("faa_path", help="Output path for .faa file")
    parser.add_argument("gff_path", help="Output path for .gff file")
    parser.add_argument("gbk_path", help="Output path for .gbk file")
    parser.add_argument("fna_path", help="Output path for .fna file")
    parser.add_argument("tsv_path", help="Output path for .tsv file")
    parser.add_argument("--meta", action="store_true", help="Metagenomic mode")
    parser.add_argument("--closed", action="store_true", help="Closed ends")
    parser.add_argument("--translation_table", type=int, default=11, help="Translation table to use")
    parser.add_argument("--write_faa", action="store_true", help="Write .faa file")
    parser.add_argument("--write_gff", action="store_true", help="Write .gff file")
    parser.add_argument("--write_gbk", action="store_true", help="Write .gbk file")
    parser.add_argument("--write_fna", action="store_true", help="Write .fna file")
    parser.add_argument("--write_tsv", action="store_true", help="Write .tsv file")
    a = parser.parse_args()

    # Run process_genome with the given arguments
    process_genome(a.sample_id, a.input_fasta, 
                 a.faa_path, a.gff_path, a.gbk_path, a.fna_path, a.tsv_path,
                 a.meta, a.closed, a.translation_table, 
                 a.write_faa, a.write_gff, a.write_gbk, a.write_fna, a.write_tsv)