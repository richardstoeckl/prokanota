# General configuration

To configure this workflow, modify `config/config.yaml` and `config/metadata.csv` according to your needs, following the explanations provided in the file.

## Pipeline configuration

To successfully run the annotation pipeline, you will need to configure the paths to the databases and the output directories in the `config/config.yaml` file.

### Optional feature prediction settings

You can optionally add a `features` section in `config/config.yaml` to control global feature prediction behavior:

```yaml
features:
	translation_table: 11      # allowed range: 1-25
	run_rrna: true             # run pybarrnap rRNA prediction
	run_trna: true             # run tRNAscan-SE prediction
	run_crispr: true           # run DICED CRISPR detection
	minimum_gene_length: 90    # minimum CDS length (bp)
```

Gene calling is selected automatically from total input length: assemblies of at least 20 kbp use single-genome training, while shorter assemblies use Pyrodigal's metagenomic models and therefore require translation table 11.
This fallback supports short assemblies; Prokanota does not currently support metagenome annotation. All contigs are treated as linear fragments, including edge CDS, until topology can be supplied per sample.
Pyrodigal masking is enabled, so runs of at least 50 unknown bases are masked during gene calling.

# Sample setup

The sample setup is specified via comma-separated tabular file (`.csv`).
Missing values can be specified by empty columns.

# Mapping file setup

Database mapping files are tab-separated (`.tsv`) with no header and exactly four columns:
`accession<TAB>short_name<TAB>description<TAB>category`.

- `accession` is the mapping key and must not be empty.
- In `short_name`, `description`, and `category`, the following values are treated as empty:
	empty/whitespace-only fields, `NA`, `N/A`, `NULL`/`null`, `-`, and `*`.
- These empty representations are normalized to `*` in final annotation outputs.

## Sample sheet

The default sample sheet is `config/metadata.csv` (as configured in `config/config.yaml`).
Each row usually corresponds to one genome with the first collumn being the sampleID (Note: The sampleID, input type, and ordered sequences are used to calculate the gene_ids!), and the second column being the path to the genome assembly in `FASTA` format. The third column indicates wether the fasta file is a genome ("dna") and the feature prediction module is run, or the fasta file is a proteome ("protein") and the feature prediction module is skipped.

## Trusted protein references

Add at most one enabled trusted entry to the databases YAML:

```yaml
databases:
  - name: trusted_protein
    search_tool: trusted_protein
    db_path: references/proteins.gbk
```

`trusted_protein` selects the specialized workflow. `db_path` is a raw GenBank or Bakta short-form protein FASTA file, not a prepared database (the parsing and building of the db is done automatically).
Relative reference paths resolve against the databases YAML directory. The entry is shared by all samples, including protein-input samples. Its name can be any valid database name. Omitted `enabled` means true; set it to false to disable the entry.
Only `name`, `search_tool`, `db_path`, and optional `enabled` are accepted.
Do not specify `order`, mapping fields, columns, thresholds, or other tool settings.
Ordinary databases can run alongside it, or it can be the only enabled "database".

FASTA headers must use Bakta's `>reference_id gene~~~product~~~dbxrefs` short-form header format. Gene and dbxrefs may be empty; product is required. Dbxrefs are ignored without validation. Duplicate reference IDs are rejected. Duplicate sequences are rejected.

GenBank import considers only CDSs without `pseudo` or `pseudogene` qualifiers. CDSs missing locus_tag, a nonempty product, or translation are skipped with a warning. Gene is optional; `hypothetical protein` is retained. There is no protein_id or translation fallback, and no dbxref/EC transfer.
If no usable proteins remain, preparation fails before feature prediction, database building, or searches.

Following Bakta's conventions, the DIAMOND environment from the DIAMOND search tool is used with `blastp --sensitive --max-target-seqs 1 --block-size 4 --index-chunks 1`. The reported best hit must meet inclusive thresholds of 90% amino-acid identity, 80% query coverage, and 80% reference coverage. Both coverages use alignment length divided by protein length.
Following Bakta’s policy, only the top DIAMOND hit is considered; if it fails the identity or coverage thresholds, no alternative hit is used.

Only per-sample `annotation/<sample>_finalAnnotation.tsv` and
`common/annotation/finalAnnotation.tsv` receive "trusted protein" columns, immediately after the base feature columns, before ordinary database columns.

- `trusted_reference_id`, `trusted_gene`, `trusted_product`

## Path resolution

Both `prokanota run` and direct Snakemake execution use these rules:

- Paths inside `config.yaml` (`global.metadata`, `global.databases`, `global.logs`, `global.interim`, `global.results`) resolve against that config file's directory.
- Database `db_path` and `mapping_path` resolve against the databases YAML directory. This includes prepared database prefixes as well as trusted raw references.
- Sample `path` values resolve against the metadata CSV directory.
- Command-line paths, including `--cfg`, `--db`, and `--meta`, resolve against the working directory. Direct Snakemake `--config` path overrides do likewise.
- Absolute paths remain absolute.
