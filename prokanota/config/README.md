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

Gene calling is selected automatically from total input length: assemblies of
at least 20 kbp use single-genome training, while shorter assemblies use
Pyrodigal's metagenomic models and therefore require translation table 11.
This fallback supports short assemblies; Prokanota does not currently support
metagenome annotation. All contigs are treated as linear fragments, including
edge CDS, until topology can be supplied per sample. Pyrodigal masking is
enabled, so runs of at least 50 unknown bases are masked during gene calling.
If omitted, these defaults are used automatically.

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
Each row usually corresponds to one genome with the first collumn being the sampleID (Note: The sampleID is used to name the output files and is used to calculate the gene_ids!), and the second column being the path to the genome assembly in `FASTA` format. The third column indicates wether the fasta file is a genome ("dna") and the feature prediction module is run, or the fasta file is a proteome ("protein") and the feature prediction module is skipped.
