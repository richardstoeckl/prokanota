# Snakemake workflow: `Prokanota`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.27.1-brightgreen.svg)](https://snakemake.github.io)

Flexible [Snakemake](https://snakemake.github.io) pipeline for **proka**ryotic **an**n**ota**tion.

## About
There are great tools available for (mostly bacterial) prokaryotic genome annotation. [Prokka](https://github.com/tseemann/prokka) has been the de facto default for many years, and [Bakta](https://github.com/oschwengers/bakta) is not only great overall, but gets regular updates making it even better.

Unfortunately, at least in my experience, annotation for archaeal species is often not perfect when done with these tools that are specialized for bacteria. Also, while the hierarchical annotation employed by these tools creates annotations that are easy to interpret, sometimes nuances get lost with them.

This is why I developed this pipeline. Inspired by the annotation approach by [Dombrowski, Nina, et al. Nature Communications, vol. 11, no. 1, Aug. 2020, p. 3939](https://doi.org/10.1038/s41467-020-17408-w), this snakemake pipeline returns the top hits from multiple databases for each predicted protein, allowing the user to make more informed decisions which annotation might be the best.

## Features

1. To ensure that the results are consistent and changes easily identifiable, the gene_ids are created by calculating a hash from the DNA content supplied genome AND the supplied sampleID. Check the "FAQ" section for more information.
2. CDS are predicted using [pyrodigal](https://github.com/althonos/pyrodigal), which fixes some bugs from prodigal. Results can differ very slightly in these edge-cases.
3. rRNAs genes are predicted using [pybarrnap](https://github.com/moshi4/pybarrnap) in `--accurate` mode. This uses the CM models from Rfam(14.10).
4. tRNAs are predicted using [tRNAscan-SE](https://github.com/UCSC-LoweLab/tRNAscan-SE) in `-G` mode.
5. To facilitate use with pangenome analyses or similar analysis, the pipeline outputs the predicted gene positions as `.gff`,`.tsv`, and `.gbk` files; the gene sequences as `.fna` files; and the protein sequences as `.faa` files. All files utilize the exact same gene_ids. Use these files as input for your favorite downstream tool and cross-reference your results with the annotations from this pipeline!
6. Currently, the genes are annotated using the [CDD](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#NCBI_curated_domains), [COG](https://www.ncbi.nlm.nih.gov/COG/), [arCOG](https://pubmed.ncbi.nlm.nih.gov/25764277/), and [PGAP](https://ftp.ncbi.nlm.nih.gov/hmm/) databases. More are planned in the future.

## Usage

**[Check out the usage instructions in the snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/richardstoeckl%20prokanota.html)**

But here is a rough overview:
1. Install [conda](https://docs.conda.io/en/latest/miniconda.html) (miniforge or miniconda is fine).
2. Install snakemake with:
```bash
conda install -c conda-forge -c bioconda snakemake
```
3. [Download the latest release from this repo](https://github.com/richardstoeckl/prokanota/releases/latest) and cd into it, or download the development version [directly from github](https://github.com/richardstoeckl/prokanota/archive/refs/heads/main.zip)
4. Edit the `config/config.yaml` to provide the paths to your results/logs directories, and the path to where you want the databases to be downloaded to.
5. Edit the `config/metadata.csv` file with the specific details for each assembly you want to annotate. Please note, that the sampleID you enter here will influence the naming of the contig and gene IDs!
5. Open a terminal in the main dir and start a dry-run of the pipeline with the following command. This will show you if you set up the paths correctly:

```bash
snakemake --sdm conda -n --cores
```
6. Run the pipeline with
```bash
snakemake --sdm conda --cores
```

## Databases and Tools used in the pipeline and reasoning. Consider citing these if you use this pipeline.
- **Databases**
    - **[CDD](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#NCBI_curated_domains) - The Conserved Domains Database curated by NCBI.** *Sources information from the [SMART](http://smart.embl-heidelberg.de/), [Pfam](http://pfam.sanger.ac.uk/), [COG](https://www.ncbi.nlm.nih.gov/COG/), [TIGR](https://www.ncbi.nlm.nih.gov/Structure/cdd/docs/tigrfams.html), and [PRK](https://www.ncbi.nlm.nih.gov/proteinclusters) databases. This makes it a good consolidated database for annotation purposes.*
    - **[COG](https://www.ncbi.nlm.nih.gov/COG/) - The Clusters of Orthologous Genes (COG) database.** *Even though the CDD sources some information from this DB as well, this is such a good resource that it deserves to be included in full and on its own.*
    - **[arCOG](https://pubmed.ncbi.nlm.nih.gov/25764277/) - The archaea specific version of the COG database.** *This is included to serve more up-to-date and accurate annotation of archaeal genomes.*
    - **[PGAP](https://ftp.ncbi.nlm.nih.gov/hmm/) - HMMs are used by the [NCBI Prokaryotic Genome Annotation Pipeline (PGAP)](https://pubmed.ncbi.nlm.nih.gov/33270901/).** *Since most submissions of novel prokaryotic genomes are evaluated by the PGAP tool, it makes sense to include its DB.*



## FAQ

- **How is the gene_id determined? What do the letters mean?**
  - The gene_id is determined by the contig_id, which is determined by hashing the DNA sequence of the genome and the user supplied sampleID. This ensures that the gene_id is consistent across runs. This is important for downstream analyses, as it allows you to cross-reference your results with the annotations from this pipeline. Changes in the contig_id indicate either changes in the sampleID or changes in the genome sequence and should probably investigated. The hashing algorithm is found in [`prokanota/workflow/scripts/cds.py`](https://github.com/richardstoeckl/prokanota/blob/97462102247b8059d4556bd9fd16e7bc85e4714d/prokanota/workflow/scripts/cds.py#L58-L84).
- **How is the predicted protein molecular weight estimated?**
  - The molecular weight is estimated by summing the average residue masses in the protein sequence and adding the mass of one water molecule, as described in [GASTEIGER, Elisabeth, et al. The proteomics protocols handbook, 2005, S. 571-607](https://doi.org/10.1385/1-59259-890-0:571). The molecular weights are taken from the [here](https://web.expasy.org/findmod/findmod_masses.html#AA) and the function can be found in [`prokanota/workflow/scripts/cds.py`](https://github.com/richardstoeckl/prokanota/blob/97462102247b8059d4556bd9fd16e7bc85e4714d/prokanota/workflow/scripts/cds.py#L240-L258)





```
Copyright Richard Stöckl 2025.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE or copy at 
https://www.boost.org/LICENSE_1_0.txt)
```