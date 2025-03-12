# General configuration

To configure this workflow, modify `config/config.yaml` and `config/metadata.csv` according to your needs, following the explanations provided in the file.

## Pipeline configuration

To successfully run the annotation pipeline, you will need to configure the paths to the databases and the output directories in the `config/config.yaml` file.

# Sample setup

The sample setup is specified via comma-separated tabular file (`.csv`).
Missing values can be specified by empty columns.

## Sample sheet

The default sample sheet is `config/metadata.csv` (as configured in `config/config.yaml`).
Each row usually corresponds to one genome with the first collumn being the sampleID (Note: The sampleID is used to name the output files and is used to calculate the gene_ids!), and the second column being the path to the genome assembly in `FASTA` format.