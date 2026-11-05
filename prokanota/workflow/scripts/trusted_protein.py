import argparse
import csv
import json
import logging
import subprocess
from pathlib import Path

from Bio import SeqIO
from io_utils import check_source_proteins, open_text, read_fasta, write_fasta_record

logger = logging.getLogger(__name__)
TRUSTED_COLUMNS = [
    "trusted_reference_id",
    "trusted_gene",
    "trusted_product",
    "trusted_identity",
    "trusted_query_coverage",
    "trusted_reference_coverage",
    "trusted_score",
    "trusted_evalue",
]


def prepare(reference, fasta, annotations):

    with open_text(reference) as handle:
        line_number, first = next(
            (
                (number, line.strip())
                for number, line in enumerate(handle, 1)
                if line.strip()
            ),
            (1, ""),
        )
    if not first or not first.startswith((">", "LOCUS")):
        raise ValueError(
            f"Reference {reference}, line {line_number}: expected nonempty Bakta short-form FASTA or GenBank"
        )
    is_fasta = first.startswith(">")
    if is_fasta:
        records = read_fasta(reference)
        check_source_proteins(reference, records)
    else:
        with open_text(reference) as handle:
            records = list(SeqIO.parse(handle, "genbank"))
    seen = set()
    sequences = {}
    mapping = {}
    output_records = []
    for record in records:
        proteins = []
        if is_fasta:
            header = record.description.split(maxsplit=1)
            fields = header[1].split("~~~") if len(header) == 2 else []
            if len(fields) != 3 or not fields[1].strip():
                raise ValueError(
                    f"Reference {reference}, record {record.id!r}, line "
                    f"{record.annotations['fasta_line']}: expected gene~~~product~~~dbxrefs "
                    "with nonempty product; extended headers are unsupported"
                )
            proteins.append(
                (record.id, fields[0].strip(), fields[1].strip(), str(record.seq))
            )
        else:
            for index, feature in enumerate(record.features, 1):
                q = feature.qualifiers
                if feature.type != "CDS" or "pseudo" in q or "pseudogene" in q:
                    continue
                accession = q.get("locus_tag", [""])[0].strip()
                product = q.get("product", [""])[0].strip()
                sequence = q.get("translation", [""])[0]
                if not accession or not product or not sequence.strip():
                    logger.warning(
                        "Skipping %s CDS %s (%s): missing locus_tag, product or translation",
                        record.id,
                        index,
                        accession,
                    )
                    continue
                proteins.append(
                    (accession, q.get("gene", [""])[0].strip(), product, sequence)
                )
        for accession, gene, product, sequence in proteins:
            if accession in seen:
                raise ValueError(
                    f"Duplicate trusted reference ID in {reference}: {accession!r}"
                )
            seen.add(accession)
            sequence = "".join(sequence.split()).upper().removesuffix("*")
            if not sequence:
                logger.warning(
                    "Skipping reference %r: empty protein sequence", accession
                )
                continue
            if sequence in sequences:
                location = (
                    f", line {record.annotations['fasta_line']}" if is_fasta else ""
                )
                raise ValueError(
                    f"Duplicate trusted protein sequence in {reference}{location}: {sequences[sequence]!r} and {accession!r}"
                )
            sequences[sequence] = accession
            # Opaque search IDs protect original accessions from DIAMOND ID parsing.
            target = f"ref{len(mapping) + 1}"
            mapping[target] = [accession, gene, product]
            output_records.append((target, sequence))
    if not mapping:
        raise ValueError(f"Trusted reference {reference} yields no usable proteins")
    with open(fasta, "w", encoding="utf-8") as output:
        for target, sequence in output_records:
            write_fasta_record(output, target, sequence)
    Path(annotations).write_text(json.dumps(mapping) + "\n", encoding="utf-8")
    logger.info("Imported %s trusted proteins from %s", len(mapping), reference)


def search(fasta, database, annotations, hits, output, version, threads):
    mapping = json.loads(Path(annotations).read_text())
    if Path(fasta).stat().st_size:
        tool_version = subprocess.check_output(
            ["diamond", "version"], text=True
        ).strip()
        Path(version).write_text(tool_version + "\n")
        subprocess.run(
            [
                "diamond",
                "blastp",
                "--query",
                fasta,
                "--db",
                database,
                "--out",
                hits,
                "--threads",
                str(threads),
                "--sensitive",
                "--load-threads",
                str(min(threads, 8)),
                "--block-size",
                "4",
                "--index-chunks",
                "1",
                "--max-target-seqs",
                "1",
                "--outfmt",
                "6",
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "qlen",
                "slen",
                "bitscore",
                "evalue",
            ],
            check=True,
        )
    else:
        Path(hits).write_text("")
        Path(version).write_text("DIAMOND skipped: empty query\n")
    with open(output, "w", newline="") as handle, open(hits) as hit_handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_id", *TRUSTED_COLUMNS])
        seen = set()
        for row in csv.reader(hit_handle, delimiter="\t"):
            query, target, identity, length, qlen, slen, score, evalue = row
            # Only the tool's first reported hit is considered, even if it fails.
            if query in seen:
                continue
            seen.add(query)
            identity, length, qlen, slen = map(float, (identity, length, qlen, slen))
            if identity >= 90 and length >= 0.8 * qlen and length >= 0.8 * slen:
                writer.writerow(
                    [
                        query,
                        *mapping[target],
                        identity,
                        100 * length / qlen,
                        100 * length / slen,
                        score,
                        evalue,
                    ]
                )


def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prep = commands.add_parser("prepare")
    for arg in ("reference", "fasta", "annotations"):
        prep.add_argument(f"--{arg}", required=True)
    run = commands.add_parser("search")
    for arg in ("fasta", "database", "annotations", "hits", "output", "version"):
        run.add_argument(f"--{arg}", required=True)
    run.add_argument("--threads", type=int, default=1)
    args = vars(parser.parse_args())
    command = args.pop("command")
    {"prepare": prepare, "search": search}[command](**args)


if __name__ == "__main__":
    main()
